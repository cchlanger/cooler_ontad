# OnTAD for coolers 
cooler_ontad is a tool that wraps [OnTAD](https://github.com/anlin00007/OnTAD) to work with the [cooler ecosystem](
https://github.com/mirnylab/cooler).
This tool consists of three steps. It creates a dense matrix for OnTAD out of a mcooler, runs OnTAD, and creates a .bedpe out of the ONTAD results, all parameters used will be extended to the filename.

## Install with pip:
pip install ${githubpath}
This will set up a command called cooler_ontad in your environment.

## Usage Examples:
This will show all possible options, the tools support all options of OnTAD:
cooler_ontad --help
This will run with the default parameters:
cooler_ontad --binsize testdata/G2.fc_1_2.cis.1000.mcool
This will run with the binsize reduced to 30000 and only add non default parameters to the filname:
cooler_ontad --shortname --binsize 30000 testdata/G2.fc_1_2.wOldG2.cis.1000.mcool

## Dockerhub container:
This tool can be found ready to use in our main container on [dockerhub](https://hub.docker.com/repository/docker/gerlichlab/scshic_docker)
