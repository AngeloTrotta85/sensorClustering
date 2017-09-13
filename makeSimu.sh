#! /bin/bash

OUTPUT_DIR="stats"
EXEC="Release/SensorsClustering"
N_RUNS=$1
MAX_OPT_OP=$2

rm -rf "$OUTPUT_DIR/algo1corr_OK_3D.data"
rm -rf "$OUTPUT_DIR/algo1corr_OK_3D.data"

rm -rf "$OUTPUT_DIR/randcorr_OK_3D.data"
rm -rf "$OUTPUT_DIR/randcorr_OK_3D.data"

rm -rf "$OUTPUT_DIR/optcorr_OK_3D.data"
rm -rf "$OUTPUT_DIR/optcorr_OK_3D.data"


#for lambda in 3 4 5
for lambda in {3..25}
do
	MIN_SENSORS=`echo $(( lambda * 2 ))`
	#echo "MIN_SENSORS $MIN_SENSORS"
	
	for sensors in {0..150..2}
	#for (( sensors=MIN_SENSORS; sensors<=120; sensors+=2 ))
	do
		if [ $sensors -lt $MIN_SENSORS ]
		then
			continue
		fi
		#echo "s:$sensors l:$lambda"
		#continue
		
		rm -rf "$OUTPUT_DIR/algo1corr_${lambda}_${sensors}.data"
		rm -rf "$OUTPUT_DIR/algo1dist_${lambda}_${sensors}.data"
		
		rm -rf "$OUTPUT_DIR/optcorr_${lambda}_${sensors}.data"
		rm -rf "$OUTPUT_DIR/optdist_${lambda}_${sensors}.data"
		
		rm -rf "$OUTPUT_DIR/randcorr_${lambda}_${sensors}.data"
		rm -rf "$OUTPUT_DIR/randdist_${lambda}_${sensors}.data"
		
		#for runs in {1..$N_RUNS}
		for (( runs=1; runs<=N_RUNS; runs++ ))
		do
			$EXEC -f inputScenario/inputTest.dat -g inputScenario/inputTest.dat -s 200 -e 0 -n $sensors -l 1 &>/dev/null
		
			NOW_T=`date +"%F %T"`
			echo "$NOW_T s:$sensors l:$lambda r:$runs"
		
			echo -n "Algo1... "
			ALGO1=`$EXEC -f inputScenario/inputTest.dat -e 4 -l $lambda | grep StatMaxCorr`
			ALGO1CORR=`echo $ALGO1 | awk '{printf $2}'`
			ALGO1DIST=`echo $ALGO1 | awk '{printf $3}'`
			echo "$ALGO1CORR" >> "$OUTPUT_DIR/algo1corr_${lambda}_${sensors}.data"
			echo "$ALGO1DIST" >> "$OUTPUT_DIR/algo1dist_${lambda}_${sensors}.data"
			#echo "Algo1 corr:$ALGO1CORR dist:$ALGO1DIST"
						
			echo -n "OK - Opt... "
			OPT=`$EXEC -f inputScenario/inputTest.dat -t $MAX_OPT_OP -e 5 -l $lambda | grep StatMaxCorr`
			OPTCORR=`echo $OPT | awk '{printf $2}'`
			OPTDIST=`echo $OPT | awk '{printf $3}'`
			echo "$OPTCORR" >> "$OUTPUT_DIR/optcorr_${lambda}_${sensors}.data"
			echo "$OPTDIST" >> "$OUTPUT_DIR/optdist_${lambda}_${sensors}.data"
			
			#EXPECTED_OP=`echo $(( (sensors/lambda) ** sensors ))`
			#if [ $EXPECTED_OP -lt $MAX_OPT_OP ]
			#then
			#	echo -n "OK - Opt ($EXPECTED_OP expected operation)... "
			#	OPT=`$EXEC -f inputScenario/inputTest.dat -e 5 -l $lambda | grep StatMaxCorr`
			#	OPTCORR=`echo $OPT | awk '{printf $2}'`
			#	OPTDIST=`echo $OPT | awk '{printf $3}'`
			#	echo "$OPTCORR" >> "$OUTPUT_DIR/optcorr_${lambda}_${sensors}.data"
			#	echo "$OPTDIST" >> "$OUTPUT_DIR/optdist_${lambda}_${sensors}.data"
			#else
			#	echo "0" >> "$OUTPUT_DIR/optcorr_${lambda}_${sensors}.data"
			#	echo "0" >> "$OUTPUT_DIR/optdist_${lambda}_${sensors}.data"
			#	echo -n "OK - Skipping Opt ($EXPECTED_OP expected operation) "
			#fi
		
			
		
			echo -n "OK - Random... "
			RAND=`$EXEC -f inputScenario/inputTest.dat -e 7 -l $lambda | grep StatMaxCorr`
			RANDCORR=`echo $RAND | awk '{printf $2}'`
			RANDDIST=`echo $RAND | awk '{printf $3}'`
			echo "$RANDCORR" >> "$OUTPUT_DIR/randcorr_${lambda}_${sensors}.data"
			echo "$RANDDIST" >> "$OUTPUT_DIR/randdist_${lambda}_${sensors}.data"
		
			echo "OK"
		done
		
		ALGO1CORR_OK=`cat "$OUTPUT_DIR/algo1corr_${lambda}_${sensors}.data" | awk 'BEGIN{c=0;s=0}{c++;s+=$1}END{print s/c}'`
		ALGO1DIST_OK=`cat "$OUTPUT_DIR/algo1dist_${lambda}_${sensors}.data" | awk 'BEGIN{c=0;s=0}{c++;s+=$1}END{print s/c}'`
		
		OPTCORR_OK=`cat "$OUTPUT_DIR/optcorr_${lambda}_${sensors}.data" | awk 'BEGIN{c=0;s=0}{c++;s+=$1}END{print s/c}'`
		OPTDIST_OK=`cat "$OUTPUT_DIR/optdist_${lambda}_${sensors}.data" | awk 'BEGIN{c=0;s=0}{c++;s+=$1}END{print s/c}'`
		
		RANDCORR_OK=`cat "$OUTPUT_DIR/randcorr_${lambda}_${sensors}.data" | awk 'BEGIN{c=0;s=0}{c++;s+=$1}END{print s/c}'`
		RANDDIST_OK=`cat "$OUTPUT_DIR/randdist_${lambda}_${sensors}.data" | awk 'BEGIN{c=0;s=0}{c++;s+=$1}END{print s/c}'`
		
		rm -rf "$OUTPUT_DIR/algo1corr_${lambda}_${sensors}.data"
		rm -rf "$OUTPUT_DIR/algo1dist_${lambda}_${sensors}.data"
		
		rm -rf "$OUTPUT_DIR/optcorr_${lambda}_${sensors}.data"
		rm -rf "$OUTPUT_DIR/optdist_${lambda}_${sensors}.data"
		
		rm -rf "$OUTPUT_DIR/randcorr_${lambda}_${sensors}.data"
		rm -rf "$OUTPUT_DIR/randdist_${lambda}_${sensors}.data"
				
		
		
		echo "$lambda $ALGO1CORR_OK" >> "$OUTPUT_DIR/algo1corr_OK_s${sensors}.data"
		echo "$sensors $ALGO1CORR_OK" >> "$OUTPUT_DIR/algo1corr_OK_l${lambda}.data"
		echo "$lambda $ALGO1DIST_OK" >> "$OUTPUT_DIR/algo1dist_OK_s${sensors}.data"
		echo "$sensors $ALGO1DIST_OK" >> "$OUTPUT_DIR/algo1dist_OK_l${lambda}.data"
		
		echo "$lambda $OPTCORR_OK" >> "$OUTPUT_DIR/optcorr_OK_s${sensors}.data"
		echo "$sensors $OPTCORR_OK" >> "$OUTPUT_DIR/optcorr_OK_l${lambda}.data"
		echo "$lambda $OPTDIST_OK" >> "$OUTPUT_DIR/optdist_OK_s${sensors}.data"
		echo "$sensors $OPTDIST_OK" >> "$OUTPUT_DIR/optdist_OK_l${lambda}.data"

		echo "$lambda $RANDCORR_OK" >> "$OUTPUT_DIR/randcorr_OK_s${sensors}.data"
		echo "$sensors $RANDCORR_OK" >> "$OUTPUT_DIR/randcorr_OK_l${lambda}.data"
		echo "$lambda $RANDDIST_OK" >> "$OUTPUT_DIR/randdist_OK_s${sensors}.data"
		echo "$sensors $RANDDIST_OK" >> "$OUTPUT_DIR/randdist_OK_l${lambda}.data"
		
		echo "$lambda $sensors $ALGO1CORR_OK" >> "$OUTPUT_DIR/algo1corr_OK_3D.data"
		echo "$lambda $sensors $ALGO1DIST_OK" >> "$OUTPUT_DIR/algo1dist_OK_3D.data"
		
		echo "$lambda $sensors $RANDCORR_OK" >> "$OUTPUT_DIR/randcorr_OK_3D.data"
		echo "$lambda $sensors $RANDDIST_OK" >> "$OUTPUT_DIR/randdist_OK_3D.data"
		
		echo "$lambda $sensors $OPTCORR_OK" >> "$OUTPUT_DIR/optcorr_OK_3D.data"
		echo "$lambda $sensors $OPTDIST_OK" >> "$OUTPUT_DIR/optdist_OK_3D.data"
	done
	
	echo "" >> "$OUTPUT_DIR/algo1corr_OK_3D.data"
	echo "" >> "$OUTPUT_DIR/algo1dist_OK_3D.data"
	
	echo "" >> "$OUTPUT_DIR/randcorr_OK_3D.data"
	echo "" >> "$OUTPUT_DIR/randdist_OK_3D.data"
	
	echo "" >> "$OUTPUT_DIR/optcorr_OK_3D.data"
	echo "" >> "$OUTPUT_DIR/optdist_OK_3D.data"
done

#Release/SensorsClustering -f inputScenario/inputTest.dat -g inputScenario/inputTest.dat -o outputRis/outputTest.dot -s 100 -e 7 -n 1200 -l 33

