help_message () {
	echo ""
	echo "Usage: run_megahit [options] -i inputdirectory -t threads -m memory -o output_dir "
	
	echo ""
	echo "Options:"
	echo ""
	echo "	-i STR          short read input directory"	
	echo "	-o STR          output directory"
	echo "	-t INT          number of threads (default=1)"
	echo "	-m INT	        amount of RAM available (default=10)"
	echo ""
   echo "";}



# default params
threads=1; mem=4; in=false; out=false

# load in params
OPTS=`getopt -o h:i:o:t:m: --long help,in,out,threads,mem`
if [ $? -ne 0 ]; then help_message; exit 1; fi

# loop through input params
while true; do
        case "$1" in
         -t) threads=$2; shift 2;;
			-i) in=$2; shift 2;;
			-m) mem=$2; shift 2;;
         -o) out=$2; shift 2;;        
         -h | --help) help_message; exit 1; shift 1;;
         --) help_message; exit 1; shift; break ;;
         *) break;;
        esac
done


   # Assembly with spades
   finallist=$(ls -d $in/*_1.fastq | awk '{print $NF}')
   
   for i in ${finallist}  
   do
      b=`basename $i _1.fastq`
     
      echo "Running megahit for sample : " $b
      mkdir -p ${out}/assembly/megahit_individal/
      cd ${out}/assembly/megahit_individal/
      megahit \
      -1 ${i} \
      -2 $in${b}_2.fastq \
      -r ${b}.fastq
      -t ${threads} \
      -o $b

   done
