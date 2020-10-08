
PREFIX=~/opt

install_prerequisites:
	# sratoolkit
	git clone git@github.com:ncbi/ngs.git
	pushd ngs
	./configure --prefix=${PREFIX}
	make
	make install
	popd

	pushd ngs/ngs-sdk
	./configure --prefix=${PREFIX}
	make
	make install
	popd

	git clone git@github.com:ncbi/ncbi-vdb.git
	pushd ncbi-vdb
	./configure --prefix=${PREFIX}
	make
	make install
	popd

	git clone git@github.com:ncbi/sra-tools.git
	cd sra-tools
	./configure --prefix=${PREFIX}
	make
	make install
	popd

	# bowtie2
	wget https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.4.2/bowtie2-2.4.2-macos-x86_64.zip
	unzip bowtie2-2.4.2-macos-x86_64.zip
	pushd bin
	ln -s ../bowtie2-2.4.2-macos-x86_64/bowtie2* .
	popd

	# RSEM
	git clone git@github.com:deweylab/RSEM.git
	pushd RSEM
	make
	make install DESTDIR=${PREFIX} prefix=${PREFIX}
	popd

make ca_coexp:
	Rscript scripts/1_define_runs.R
	Rscript scripts/2_download_runs.R
	Rscript scripts/3_reference_genome.R
	Rscript scripts/4_submit_estimate_expression.R
