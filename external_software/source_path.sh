#### Get the directory of the script:
export EXTERNAL_SOFTWARE_DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

echo -e "\n\n"
type conda >/dev/null 2>&1 || { echo >&2 "Conda is required to be installed... and is not being found in the PATH. I'm assuming conda has been installed within reanalyzerGSE/external_software folder and adding to the path. Otherwise, please double check manually or install conda beforehand..."; export conda_within_reanalyzer="yes"; }

if [ "$conda_within_reanalyzer" == "yes" ]; then
	conda_dir=$EXTERNAL_SOFTWARE_DIR/miniconda3/bin
else	
	conda_dir=$(dirname $(which conda))
fi

conda_envs_path=$(echo $conda_dir | sed 's,/bin$,/envs,g')

#### Set PATH:

echo -e "\nPlease run manually next and before reanalyzerGSE execution:\nexport PATH=$conda_envs_path/reanalyzerGSE/bin:$conda_dir:$EXTERNAL_SOFTWARE_DIR/miARma-seq:$(dirname $EXTERNAL_SOFTWARE_DIR):$(dirname $EXTERNAL_SOFTWARE_DIR)/scripts:$HOME/bin:$PATH\n"

echo -e "\n\nPlease behave and be mindful with any queueing system, the resources you are using, the versions of the software that are in the PATH and are being used... etc\n\n"