# variant_liftover
A basic tool for lit over variant call files from one genome to another without need of liftover chain

# Requirements
The tool requires samtools and bowtie2. Both tools are part of a docker container used in our metagenome pipeline
that can be recycled here. To run on a hpc, e.g run

```
singularity build assemble03.sif docker://fischuu/hrp_assemble:0.3
singularity shell assemble03.sif
```

And then run the tool in that container.

