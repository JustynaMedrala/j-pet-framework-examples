filename="/mnt/home/jmedrala/output/mcGeant/1_larger/1695_2024_07_01-04_04_23.mcGeant.root"
base=$(basename ${filename})
echo "${{$(basename ${filename})}%.*.*}.hits.evt.root"
