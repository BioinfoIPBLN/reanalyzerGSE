#!/bin/bash
echo -e "\nThis is a wrapper script to install all reanalyzerGSE dependencies available through conda (multiple environments due to dependency conflicts), and others that are external, such as check_strandedness with pip.\n"
echo -e "\nBy default, the script first checks if there's a conda executable in the PATH. If there is, it tries to install mamba in the active environment (it shouldn't change anything if already installed) and then reads four yml files to create new environments. If there is no conda installed, it tries to install latest conda version in the external_software folder and performs the same steps in the sentence before\n"
echo -e "\nIf you want to use just conda instead of mamba, please use '-m conda'\n"

### Define arguments and variables
### A string with command options and an array with arguments
options=$@
arguments=($options)
### Get arguments by looping through an index
index=0
for argument in $options; do
### Incrementing index
	index=`expr $index + 1`
### Gather the parameters
	case $argument in
		-h*) echo "install.sh usage: install.sh [options]
		-h | -help # Type this to get help
		-m | -manager # manager to use (mamba, by default, or conda)" && exit 1;;
		-m) manager=${arguments[index]}
	esac
done

if [ -z "$manager" ]; then
	manager="mamba"
fi


#### Get the directory of the script:
CURRENT_DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
echo -e "\nThe current folder is $CURRENT_DIR\n"

#### Setting permissions of executables:
echo -e "\n\nSetting permissions...\n\n"
chmod 775 $CURRENT_DIR/../../reanalyzerGSE.pk.sh
chmod 775 $CURRENT_DIR/../../scripts/*
chmod 775 $CURRENT_DIR/../miARma-seq/miARma

#### Download and install Miniconda if required:
echo -e "\n\nInstalling dependencies through conda...\n\n"
export conda_install=""
type conda >/dev/null 2>&1 || { echo >&2 "Conda is required to be installed... and is not being found in the PATH. Conda installation within reanalyzerGSE/external_software folder is going to happen in 30 seconds. Please kill the process and correct the PATH if not necessary, otherwise let the installation continue..."; export conda_install="yes"; }
secs=$((1 * 30))
while [ $secs -gt 0 ]; do
   echo -ne "$secs\033[0K\r"
   sleep 1
   : $((secs--))
done

if [ "$conda_install" == "yes" ]; then
	echo -e "\n\n\nI'm downloading and installing an updated copy of Miniconda3 in the folder external_software/miniconda3 to create environments...\n\n\n"
	echo -e "\nIf for some reason an outdated python or conda is required in your system, please kill this process and and go to https://repo.anaconda.com/miniconda/ to download and install manually the corresponding Linux installer\n"
	echo -e "\nProceeding with conda installation\n"
	cd $CURRENT_DIR/.. && wget -q "https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh"
	mkdir -p miniconda3 && bash Miniconda3-latest-Linux-x86_64.sh -b -f -s -p $PWD/miniconda3 && rm Miniconda3-latest-Linux-x86_64.sh
	export PATH=$CURRENT_DIR/../miniconda3/bin:$PATH
	if [[ ! -f $CURRENT_DIR/../miniconda3/bin/conda ]]; then
		echo "Installation of miniconda3 failed, please check manually or install manually and put conda in the PATH... Exiting..."
   		exit 1
	fi
fi

#### Install packages via conda:   
conda_exec=$(which conda)
conda_dir=$(dirname $conda_exec | sed 's,/bin$,,g' | sed 's,/condabin$,,g')
conda_envs_path=$(dirname $conda_exec | sed 's,/bin$,/envs,g' | sed 's,/condabin$,/envs,g')
echo -e "\n\n\nI'm downloading and installing several packages in four environments in $conda_dir...\nThe conda executable used is located in the conda located in $conda_exec\n\n"
echo -e "\n\n\nThe pathway to conda environments is $conda_envs_path...\n\n\n"

if [[ $manager == "mamba" ]]; then
	echo -e "The manager chosen is 'mamba', first trying to install mamba in the base environment replacing 'conda' to reduce time, if it's not already installed, and then populating the environments based on the yml files (keep in mind, this is frozen versions of the software)...\n\n\n"
	conda install -y -q -c conda-forge mamba
 	echo -e "\n\nInstalling reanalyzerGSE1...\n\n"
	mamba env create -q --file $CURRENT_DIR/reanalyzerGSE1.yml
 	echo -e "\n\nInstalling reanalyzerGSE2...\n\n"
	mamba env create -q --file $CURRENT_DIR/reanalyzerGSE2.yml
 	echo -e "\n\nInstalling reanalyzerGSE3...\n\n"
	mamba env create -q --file $CURRENT_DIR/reanalyzerGSE3.yml
 	echo -e "\n\nInstalling reanalyzerGSE4...\n\n"
	mamba env create -q --file $CURRENT_DIR/reanalyzerGSE4.yml
else
	echo -e "The manager chosen is 'conda', populating the environments based on the yml files (keep in mind, this is frozen versions of the software)...\n\n\n"
 	echo -e "\n\nInstalling reanalyzerGSE1...\n\n"
	conda env create -q --file $CURRENT_DIR/reanalyzerGSE1.yml
 	echo -e "\n\nInstalling reanalyzerGSE2...\n\n"
	conda env create -q --file $CURRENT_DIR/reanalyzerGSE2.yml
 	echo -e "\n\nInstalling reanalyzerGSE3...\n\n"
	conda env create -q --file $CURRENT_DIR/reanalyzerGSE3.yml
 	echo -e "\n\nInstalling reanalyzerGSE4...\n\n"
	conda env create -q --file $CURRENT_DIR/reanalyzerGSE4.yml
fi

echo -e "\nRemoving tmp files in pkgs directory\n"
rm -rf $(find $(which conda | sed 's,/bin/conda,,g') -type d -name pkgs) # Remove temp files

$conda_envs_path/reanalyzerGSE_4/bin/pip -q install statistics
$conda_envs_path/reanalyzerGSE/bin/pip -q install deeptools

echo -e "\n\nSoft linking some software so only one environment/path has to be used...\n\n"
cd $conda_envs_path/reanalyzerGSE/bin
ln -sf $conda_envs_path/reanalyzerGSE_2/bin/featureCounts .
ln -sf $conda_envs_path/reanalyzerGSE_2/bin/salmon .
ln -sf $conda_envs_path/reanalyzerGSE_2/bin/seqtk .
ln -sf $conda_envs_path/reanalyzerGSE_2/bin/hisat2 .
ln -sf $conda_envs_path/reanalyzerGSE_2/bin/hisat2-build .
ln -sf $conda_envs_path/reanalyzerGSE_2/bin/taxonkit .
ln -sf $conda_envs_path/reanalyzerGSE_2/bin/extract_kraken_reads.py .
ln -sf $conda_envs_path/reanalyzerGSE_2/bin/kraken2 .
ln -sf $conda_envs_path/reanalyzerGSE_2/bin/sortmerna .
ln -sf $conda_envs_path/reanalyzerGSE_3/bin/R .
ln -sf $conda_envs_path/reanalyzerGSE_3/bin/Rscript .
ln -sf $conda_envs_path/reanalyzerGSE_3/bin/fastq-dl .
ln -sf $conda_envs_path/reanalyzerGSE_3/bin/qualimap .
ln -sf $conda_envs_path/reanalyzerGSE_3/bin/multiqc .
ln -sf $conda_envs_path/reanalyzerGSE_3/bin/pandoc .
ln -sf $conda_envs_path/reanalyzerGSE_3/bin/rcf .
ln -sf $conda_envs_path/reanalyzerGSE_4/bin/check_strandedness .
ln -sf $conda_envs_path/reanalyzerGSE_4/bin/gff32gtf .
ln -sf $conda_envs_path/reanalyzerGSE_4/bin/gtf2bed .
ln -sf $conda_envs_path/reanalyzerGSE_4/bin/infer_experiment.py .
ln -sf $conda_envs_path/reanalyzerGSE_4/bin/kallisto .

echo -e "\nSmall fix on check_strandedness...\n"
cd $(dirname $(find $conda_envs_path/reanalyzerGSE_4 -name check_strandedness.py))
rm check_strandedness.py; wget -q https://github.com/signalbash/how_are_we_stranded_here/raw/master/how_are_we_stranded_here/check_strandedness.py; chmod 775 check_strandedness.py

echo -e "\n\nInstalling the perl modules required for miARma-Seq via cpanm...\n\n"
cd $conda_envs_path/reanalyzerGSE/bin
curl -L https://cpanmin.us | $conda_envs_path/reanalyzerGSE/bin/perl - App::cpanminus
$conda_envs_path/reanalyzerGSE/bin/cpanm -q -f Config::IniFiles DateTime Strict::Perl LWP Cwd less Statistics::R

echo -e "\n\nInstalling the CRAN's package 'autoGO' and some other dependencies, not available through conda as of yet... This manual installation implies that the version of the installed R packages are not frozen. So, not likely, but please do keep in mind that this may be a source of errors in the mid-term if newer versions of the R packages are installed...\n\n"
export PATH=$conda_envs_path/reanalyzerGSE_3/bin:$PATH
$conda_envs_path/reanalyzerGSE_3/bin/Rscript -e 'suppressMessages({install.packages(c("autoGO","GOxploreR","aPEAR"),repos="https://cloud.r-project.org",quiet=T)})' &> /dev/null # Avoid the very extensive logs being printed out, please remove suppressMessages and "&> /dev/null" to double check if installation fails...
$conda_envs_path/reanalyzerGSE_3/bin/Rscript -e 'packages <- c("autoGO", "GOxploreR", "aPEAR"); for (pkg in packages) { if(requireNamespace(pkg, quietly=TRUE)){cat("\n",pkg, "is installed.\n")} else {cat("\n",pkg, "is NOT installed.\n")}}'

# https://github.com/DerrickWood/kraken2/issues/518, you just have to manually replace 'ftp' by 'https' in the line 46 of the file 'rsync_from_ncbi.pl'
# Fix and improve kraken2 2.1.2 within conda (i.e. fix download for database building and improve masking with parallel):
echo -e "\n\n\nI'm fixing some issues with kraken2...(see the content and comments within the script)\n\n\n"
sed -i 's,#^ftp:,#^https:,g' $conda_envs_path/reanalyzerGSE_2/libexec/rsync_from_ncbi.pl
# Fixing that it uses the perl of the corresponding environment:
sed -i "s,/usr/bin/env perl,$conda_envs_path/reanalyzerGSE_2/bin/perl,g" $conda_envs_path/reanalyzerGSE_2/bin/kraken2
# Database building is very slow, I incorporate this suggestion to paralelize, improving with multithreading and pipepart: https://github.com/DerrickWood/kraken2/pull/39
cd $conda_envs_path/reanalyzerGSE_2/libexec/
set +H
echo -e "#!/bin/bash\n" > temp
echo "function mask_data_chunk () { MASKER=\$1; awk -v RS=\">\" -v FS=\"\n\" -v ORS=\"\" ' { if (\$2) print \">\"\$0 } ' | \$MASKER -in - -outfmt fasta | sed -e '/^>/!s/[a-z]/x/g' }; export -f mask_data_chunk" | cat - mask_low_complexity.sh |  sed 's,#!/bin/bash,,g' >> temp
sed -i '/$MASKER -in $file -outfmt fasta/c\      parallel --pipepart -j $KRAKEN2_THREAD_CT -a $file --recstart ">" --block 250M mask_data_chunk $MASKER > $file.tmp' temp
sed -i '/$MASKER -in $target -outfmt fasta/c\    parallel --pipepart -j $KRAKEN2_THREAD_CT -a $target --recstart ">" --block 250M mask_data_chunk $MASKER > $target.tmp' temp
sed -i 's,) { ,) { \n  ,g' temp
sed -i 's,; awk,\n  awk,g' temp
sed -i "s,x/g' },x/g'\n},g" temp
sed -i 's,; export,\nexport,g' temp
rm mask_low_complexity.sh; mv temp mask_low_complexity.sh; chmod 775 mask_low_complexity.sh

# Fix and improve krakentools: (downloading a new version of the script not in the release, and changing the shebang so it uses the corresponding python version)
cd $conda_envs_path/reanalyzerGSE_2/bin
rm extract_kraken_reads.py; wget -q https://github.com/jenniferlu717/KrakenTools/raw/master/extract_kraken_reads.py
# Fixing that it uses the perl of the corresponding environment:
sed -i "s,/usr/bin/env python,$conda_envs_path/reanalyzerGSE_2/bin/python,g" extract_kraken_reads.py


echo -e "\n\n\nALL DONE. Please remember that reanalyzerGSE needs the conda environment 'reanalyzerGSE' to be activated, or the path '$conda_envs_path/reanalyzerGSE/bin' added to the exported PATH variable. The script located in 'external_software/source_path.sh' will also suggest the proper PATH to be exported\n\n\n"






