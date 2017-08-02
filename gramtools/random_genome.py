import random


def _generate_random_genome(genome_regions):
    """Generate a genome by selecting a single allele from each block."""
    genome = [random.choice(region.alleles) for region in genome_regions]
    return ''.join(genome)


def dump_random_genome(nreads, genome_regions):
    random_genome = _generate_random_genome(genome_regions)

    with open("random_genome.fa", "w") as fhandle:
        fhandle.write(">randGenome\n{}\n".format(random_genome))

    fhandle = open("random_reads.fastq", "w")
    for i in range(nreads):
        posic = random.randint(0, len(random_genome) - 150)
        fhandle.write("@{0}\n{1}\n+\n{2}\n".format(
            i, random_genome[posic: posic + 150], 'H' * 150))
    fhandle.close()
