while read line;
do
intervene upset -i ${line}_AN1.txt \
		  ${line}_AN2.txt \
		  ${line}_COP.txt \
		  ${line}_CUL.txt \
		  ${line}_F.txt \
		  ${line}_I.txt \
		  ${line}_II.txt \
		  ${line}_III.txt \
		  ${line}_IX.txt \
		  ${line}_PF.txt \
		  ${line}_PRM.txt \
		  ${line}_SIM.txt \
		  ${line}_VI.txt \
		  ${line}_VII.txt \
		  ${line}_VIII.txt \
		  ${line}_X.txt \
        --type list -o out_upset_${line} --save-overlaps --figtype "png"
done < sample_ID.txt
