# Make sure you create a directory where you want to keep your log files in this case (~/log)
script_dir=$(dirname $0)
cd $script_dir/log && mv score.log score`date +"%Y_%m_%d"`.log && touch score.log