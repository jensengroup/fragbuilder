/home/andersx/programs/dsspcmbi -na tmp.pdb | grep -Ev "\.$" | grep -v "  #" | cut -b 17 | tr -d "\n" | sed -e "s/[^HE]/C/g" > tmp.sec
