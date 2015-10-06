## Transfuse

Transfuse intelligently merges your multiple de novo transcriptome assemblies. Run multiple assemblies with different de novo assemblers, or different settings in the same assembler and have them combined into a single high quality transcriptome.

Transfuse takes in the reads you used to do the assembly and a list of fasta files and produces a single output fasta file.

### Installation and Running

To install Transfuse you can get it from rubygems.com

`gem install transfuse`

or you can clone this repo:

`git clone https://github.com/cboursnell/transfuse.git`

then build and install the ruby gem

`gem build *spec; gem install *gem`

Transfuse also requires `vsearch` to be installed which can be downloaded from:

`https://github.com/torognes/vsearch`

### Usage

Transfuse is run on the command line. The options are:

```
  -a, --assemblies=<s>    assembly files in FASTA format, comma-separated
  -l, --left=<s>          left reads file in FASTQ format
  -r, --right=<s>         right reads file in FASTQ format
  -o, --output=<s>        write merged assembly to file
  -t, --threads=<i>       number of threads (default: 1)
  -i, --id=<f>            sequence identity to cluster at (default: 1.0)
  -v, --verbose           be verbose
  -e, --version           Print version and exit
  -h, --help              Show this message
```

An example command:

`transfuse --assembly soap-k31.fa,soap-k41.fa,soap-k51.fa --left reads_1.fq --right reads_2.fq --output soap-merged.fa --threads 12`

### Contributing

[![Join the chat at https://gitter.im/cboursnell/transfuse](https://badges.gitter.im/Join%20Chat.svg)](https://gitter.im/cboursnell/transfuse?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)

Tranfuse is currently in development.

If you want to suggest, and maybe implement, a new feature, please suggest it on the tracker first.

### License

This is adademic software - please cite us if you use it in your work.

Transfuse is released under the MIT license.
