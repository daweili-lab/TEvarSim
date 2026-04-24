import re
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq

def SeqDiverse(seq: str,
               snp_rate: float = 0.02,
               indel_rate: float = 0.005,
               ins_ratio: float = 0.5,
               indel_geom_p: float = 0.7,
               max_indel_len: int = 20) -> str:
    seq_array = list(seq)
    L = len(seq_array)
    BASES = ['A', 'C', 'G', 'T']
    CHOICES_DICT = {b: [c for c in BASES if c != b] for b in BASES}
    # INDEL
    n_indel = np.random.poisson(indel_rate * L)
    #print(n_indel)
    n_indel = min(n_indel, int(L * 0.1))
    #print(n_indel)
    if n_indel > 0:
        # generate INDEL lengths
        indel_len_list = np.random.geometric(indel_geom_p, n_indel)
        indel_len_list = np.clip(indel_len_list, a_min=None, a_max=int(max_indel_len))
    # select positions for INDELs
        indel_positions = np.random.choice(L, size=min(n_indel, L), replace=False)
        for pos, indel_len in zip(sorted(indel_positions, reverse=True), indel_len_list):
            if np.random.random() < ins_ratio:
                ins_seq = [np.random.choice(BASES) for _ in range(indel_len)]
                seq_array[pos:pos] = ins_seq
            else:
                del_end = min(pos + indel_len, len(seq_array))
                del seq_array[pos:del_end]
    # SNP
    L = len(seq_array)
    n_snp = np.random.poisson(snp_rate * L)
    if n_snp > 0:
        snp_positions = np.random.choice(L, size=min(n_snp, L), replace=False)
        for pos in snp_positions:
            current_base = seq_array[pos]
            if current_base in CHOICES_DICT:
                seq_array[pos] = np.random.choice(CHOICES_DICT[current_base])
    return ''.join(seq_array)

def Get_config(config_file):
    divConfig = {}
    with open(config_file, 'r') as f:
        for line in f:
            if line.strip() and not line.startswith('#'):
                key, value = line.strip().split('=')
                divConfig[key.strip()] = float(value.strip())
    return divConfig

class Simulator:
    def __init__(self, args):
        self.reference = args.ref
        self.pool_fasta = args.pool
        self.bed_file = args.bed
        self.output_prefix = args.outprefix
        self.num_genomes = args.num
        self.af_min = args.af_min
        self.af_max = args.af_max
        self.tsd_min = args.tsd_min
        self.tsd_max = args.tsd_max
        self.sense_strand_ratio = args.sense_strand_ratio
        self.diverse = args.diverse
        self.diverse_config = args.diverse_config
        self.random_seed = args.seed

        self.TEevents = []
        #self.genotypes = ""
        if self.random_seed is not None:
            np.random.seed(self.random_seed)
    
    def _run(self):
        self._parse_bed()
        self._check_bed()
        self._random_sample_genotypes()
        self.get_TE_tag()
        self.generate_vcf()
        self.generate_genome()

    def _check_bed(self):
        ref_record = next(SeqIO.parse(self.reference, "fasta"))
        self.ref_id = ref_record.id
        self.ref_len = len(ref_record.seq)
        self.CHR = ref_record.id
        self.ref_seq = str(ref_record.seq)
        prev_start = -1
        prev_end = -1
        for event in self.TEevents:
            start = event["start"]
            end = event["end"]
            if start >= end:
                raise ValueError(f"Invalid TE event: start >= end. Position: {start}")
            if start < 0 or end > len(self.ref_seq):
                raise ValueError(f"TE event out of genome bounds. Position: {start}")
            if start < prev_start:
                raise ValueError(f"bed file not be sorted by start position. Position: {start}")
            if start < prev_end:
                raise ValueError(f"Overlapping TE events detected. Position: {start}")

    def _parse_bed(self):
        """
        """
        self.TEevents = []
        with open(self.bed_file) as f:
            for line in f:
                if line.startswith("#") or not line.strip():
                    continue
                _, start, end, te_id = line.strip().split("\t")
                te_id = "DEL" if te_id == "-" else te_id # If TE deletion has no ID, assign "DEL"
                start = int(start)
                end = int(end)
                event_type = "INS" if end == start + 1 else "DEL"
                self.TEevents.append({
                    "start": start,
                    "end": end,
                    "te_id": te_id,
                    "type": event_type
                })
        print(f"[INFO] Parsed {len(self.TEevents)} TE events from BED.")
        print(f"[INFO] Example event: {self.TEevents[0] if self.TEevents else 'No events'}")
    
    def _random_sample_genotypes(self):
        """
        Randomly generate an allele frequency for each TE event, and then generate the genotype for each sample based on that frequency.
        """
        # 1. allele frequency
        if self.af_min > self.af_max:
            raise ValueError(
        f"Parameter error: af_min must be less than or equal to af_max.")
        nTE = len(self.TEevents)
        afs = np.random.uniform(self.af_min, self.af_max, size=nTE)

        # 2. sample genotypes
        self.genotypes = np.zeros((nTE, self.num_genomes), dtype=int)
        for i, af in enumerate(afs):
            self.genotypes[i] = np.random.binomial(1, af, size=self.num_genomes)

        print(f"[INFO] Generated genotypes for {nTE} events across {self.num_genomes} genomes.")
        print(f"[INFO] Allele frequencies (first 10): {afs[:10]}")
        # print(f"[INFO] Example genotypes (first 5 events):\n{self.genotypes[:5]}")

        return self.genotypes    

    def get_TE_tag(self):
        """
        Generate additional metadata for each TE event:
        strand: the insertion orientation ("+" or "-")
        tsd_len: TSD length (from tsd_min to tsd_max)
        ref: reference sequence
        alt: variant sequence (including the result after TE insertion or deletion)"
        """
        # TE pool
        te_pool = {rec.id: rec.seq for rec in SeqIO.parse(self.pool_fasta, "fasta")}
        for event in self.TEevents:
            start, end, te_id = event["start"], event["end"], event["te_id"]
            # strand
            strand = "+" if np.random.rand() < self.sense_strand_ratio else "-"
            # tsd length
            tsd_len = np.random.randint(self.tsd_min, self.tsd_max + 1)
            if te_id.startswith("bg"):
                tsd_len = 0
            if event["type"] == "INS":
                # bed file is 0-based
                ref_allele = self.ref_seq[start - 1]
                te_seq = te_pool[te_id]
                if strand == "-":
                    te_seq = str(te_seq.reverse_complement())
                # tsd_seq = ref_seq[start-tsd_len : start]
                tsd_seq = self.ref_seq[start-tsd_len+1 : start+1]
                alt_allele = self.ref_seq[start - 1] + str(te_seq) + tsd_seq
            else:
                # REF sequence
                ref_allele = self.ref_seq[start-1:end]
                alt_allele = self.ref_seq[start-1]
            # update TE tags
            event.update({
                "strand": strand,
                "tsd_len": tsd_len,
                "ref": ref_allele,
                "alt": alt_allele
            })
    
    def _parse_te_modification(self, te_id: str):
        """
        Parse TE ID into family name and modifications.
        Example: "LINE_1SNP0INDEL5polyA" → ("LINE", {"nSNP":1, "npolyA":5})
        """
        te_family = te_id
        mods = {}
        if "_" in te_id:
            te_family, mod_str = te_id.rsplit("_", 1)
            for key in ["SNP", "INDEL", "polyA", "truncate"]:
                match = re.search(rf"(\d+){key}", mod_str)
                if match:
                    val = int(match.group(1))
                    if val > 0:
                        mods["n" + key] = val
        return te_family, mods
    

    def generate_vcf(self):
        """
        Generate a VCF file from parsed TE events with detailed INFO.
        INFO includes:
            - TYPE: INS or DEL
            - STRAND: +/-
            - TSD: Target site duplication length
            - TEFAMILY: TE family ID (from consensus)
            - SNP/INDEL/polyA/truncate: modification counts (only for INS)
        """
        vcf_path = f"{self.output_prefix}.vcf"
        with open(vcf_path, "w") as vcf:
            # VCF Header
            vcf.write("##fileformat=VCFv4.2\n")
            vcf.write(f"##contig=<ID={self.ref_id},length={self.ref_len}>\n")
            vcf.write('##INFO=<ID=TYPE,Number=1,Type=String,Description="Variant type (INS or DEL)">\n')
            vcf.write('##INFO=<ID=STRAND,Number=1,Type=String,Description="Insertion strand (+ or -)">\n')
            vcf.write('##INFO=<ID=TSD,Number=1,Type=Integer,Description="Target site duplication length">\n')
            vcf.write('##INFO=<ID=TEFAMILY,Number=1,Type=String,Description="TE family ID from consensus">\n')
            vcf.write('##INFO=<ID=SNP,Number=1,Type=Integer,Description="Number of SNP modifications in TE sequence">\n')
            vcf.write('##INFO=<ID=INDEL,Number=1,Type=Integer,Description="Number of INDEL modifications in TE sequence">\n')
            vcf.write('##INFO=<ID=polyA,Number=1,Type=Integer,Description="Length of added polyA tail in TE sequence">\n')
            vcf.write('##INFO=<ID=truncate,Number=1,Type=Integer,Description="Number of truncated bases at 5\' end of TE sequence">\n')
            vcf.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
            samples = [f"Hap{i+1}" for i in range(self.num_genomes)]
            vcf.write(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples) + "\n")

            for idx, event in enumerate(self.TEevents):
                chrom = self.CHR
                # Although BED files are 0-based and VCF files are 1-based, 
                # VCF files require retrieving one base upstream, so the overall position remains unchanged.
                pos = event["start"]
                var_id = event["te_id"]
                ref = event["ref"]
                alt = str(event["alt"])
                event_type = event["type"]
                strand = event["strand"]
                tsd = event["tsd_len"]

                # parse TEID: TE family + modified information
                te_family, mods = self._parse_te_modification(var_id)

                # INFO
                info_parts = [
                    f"TYPE={event_type}",
                    f"STRAND={strand}",
                    f"TSD={tsd}",
                    f"TEFAMILY={te_family}"
                ]
                if event_type == "INS":  # only insert needs to be changed
                    for key in ["SNP", "INDEL", "polyA", "truncate"]:
                        val = mods.get("n" + key)
                        if val is not None:
                            info_parts.append(f"{key}={val}")

                info_str = ";".join(info_parts)
                genotypes = list(map(str, self.genotypes[idx]))
                vcf.write(f"{chrom}\t{pos}\t{var_id}\t{ref}\t{alt}\t.\tPASS\t{info_str}\tGT\t" + "\t".join(genotypes) + "\n")

        print(f"[INFO] VCF file written to {vcf_path}")


    def generate_genome(self):
        """
        Generate the genome sequence for each sample based on the genotype matrix and variant information.
        """
        # TE sequences
        ref_seq = str(next(SeqIO.parse(self.reference, "fasta")).seq)

        # split genome
        chunks = []
        cols_to_replace = []
        flipped = []
        col_index = 0
        start = 0
        for event in self.TEevents:
            SVtype = event['type']
            left = event['start']
            right = event['end']
            if left != start:
                if SVtype == "INS":
                    chunks.append(ref_seq[start:left+1])
                else:
                    chunks.append(ref_seq[start:left])
                col_index += 1
            if SVtype == "INS":
                chunks.append(event["alt"][1:])
            else:
                chunks.append(event["ref"][1:])
                flipped.append(col_index)
            cols_to_replace.append(col_index)
            col_index += 1
            start = right    

        # Check if the last TE occurs at the end of the genome
        genome_len = len(ref_seq)
        if right != genome_len:
            chunks.append(ref_seq[right:genome_len])

        # combine genome
        if self.diverse:
            print("introduce sequence diversity for each TE-events")

        with open(f"{self.output_prefix}.fa", "w") as fo:
            indexMat = np.ones((len(chunks), self.num_genomes))
            indexMat[cols_to_replace,:] = self.genotypes
            for idx in range(self.num_genomes):
                mask = indexMat[:, idx].astype(bool)
                mask[flipped] = ~mask[flipped]
                if self.diverse:
                    # introduce sequence diversity for each TE-events
                    mask_seq = [m if i in cols_to_replace else 0 for i, m in enumerate(mask)]
                    if self.diverse_config:
                        divConfig = Get_config(self.diverse_config)
                        diverse_chunks = [SeqDiverse(chunk, **divConfig) if use else chunk for chunk, use in zip(chunks, mask_seq)]
                    else:
                        diverse_chunks = [SeqDiverse(chunk) if use else chunk for chunk, use in zip(chunks, mask_seq)]
                    assembly = ''.join(chunk for chunk, use in zip(diverse_chunks, mask) if use)
                else:
                    assembly = ''.join(chunk for chunk, use in zip(chunks, mask) if use)
                fo.write(f">{self.ref_id}_{idx}\n")
                for i in range(0, len(assembly), 60):
                    fo.write(assembly[i:i+60] + "\n") 

        """
        for k, v in locals().items():
            if k != 'self':
                setattr(self, k, v)          
        """

def run(args):
    Simulator(args)._run()
