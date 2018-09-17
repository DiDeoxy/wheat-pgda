
# set the size of a mega base
mb <- 1000000

# open the gds object specified in the scrips
wheat <- snpgdsOpen(gds)

# extract the marker metadata
snp_id <- as.character(read.gdsn(index.gdsn(wheat, "snp.id")))
snp_chrom <- as.integer(read.gdsn(index.gdsn(wheat, "snp.chromosome")))
# turn snp base positons into Mb positons
snp_pos <- as.integer(read.gdsn(index.gdsn(wheat, "snp.position"))) / mb
# create a data frame of the marker metadata
data <- tibble(id = snp_id, chrom = snp_chrom, pos = snp_pos)

## find the max position of any marker on each genome for xlims
chrom_lengths <- by(snp_pos, snp_chrom, max)
max_genome_lengths <- c(max(chrom_lengths[seq(1, 19, 3)]), # A genome
                        max(chrom_lengths[seq(2, 20, 3)]), # B genome
                        max(chrom_lengths[seq(3, 21, 3)])) # D genome

# extract the genotypes and the samples
genotypes <- read.gdsn(index.gdsn(wheat, "genotype"))
sample_id <- as.character(read.gdsn(index.gdsn(wheat, "sample.id")))

# extract the categorical information
samp_annot <- read.gdsn(index.gdsn(wheat, "samp_annot"))
bp <- samp_annot$BP
mc <- samp_annot$MC
desig <- samp_annot$designation
texture <- samp_annot$texture
colour <- samp_annot$colour
habit <- samp_annot$habit

# transform year into binned eras
year <- samp_annot$Year
era <- cut(as.numeric(as.character(year)),
    breaks = c(1800, 1920, 1940, 1960, 1980, 2000, 2020),
    labels = c(
        "Pre-1920", "1921-1940", "1941-1960", "1961-1980", "1981-2000",
        "2001-2016"
    )
)
era <- addNA(era)
levels(era)[7] <- "UNKNOWN"

snpgdsClose(wheat)
