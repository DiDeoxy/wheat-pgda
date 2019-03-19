# convert ncbi identifiers to gene names
# resi
sed -i -e 's/query,/query,accession,/g' data/intermediate/blast/top_resi.csv
sed -i -e 's/AY325736/PM3-B,AY325736/g' data/intermediate/blast/top_resi.csv
sed -i -e 's/KF572030/PM8,KF572030/g' data/intermediate/blast/top_resi.csv
sed -i -e 's/MH077963/PM17,MH077963/g' data/intermediate/blast/top_resi.csv
sed -i -e 's/FJ876280/LR21,FJ876280/g' data/intermediate/blast/top_resi.csv
sed -i -e 's/AY270159/RGA2,AY270159/g' data/intermediate/blast/top_resi.csv #~LR10
sed -i -e 's/KF031303/SR33,KF031303/g' data/intermediate/blast/top_resi.csv
sed -i -e 's/LN883757/SR45,LN883757/g' data/intermediate/blast/top_resi.csv
sed -i -e 's/HM133634/UGT-S3,HM133634/g' data/intermediate/blast/top_resi.csv
sed -i -e 's/AK455349/UGT-CS,AK455349/g' data/intermediate/blast/top_resi.csv
sed -i -e 's/KU562861/WRKY70,KU562861/g' data/intermediate/blast/top_resi.csv
sed -i -e 's/KT962210/ACT,KT962210/g' data/intermediate/blast/top_resi.csv
sed -i -e 's/KC244204/GLI1,KC244204/g' data/intermediate/blast/top_resi.csv
sed -i -e 's/JN631792/YR5,JN631792/g' data/intermediate/blast/top_resi.csv
sed -i -e 's/AJ438338/UGT-1i15e,AJ438338/g' data/intermediate/blast/top_resi.csv
sed -i -e 's/AJ438331/UGT-3i11b,AJ438331/g' data/intermediate/blast/top_resi.csv
sed -i -e 's/FJ535237/UGT-WS-2,FJ535237/g' data/intermediate/blast/top_resi.csv
sed -i -e 's/KX274091/UGT-UNKNOWN,KX274091/g' data/intermediate/blast/top_resi.csv
sed -i -e 's/EU568801/UGT1,EU568801/g' data/intermediate/blast/top_resi.csv
sed -i -e 's/AJ438332/UGT-3i11e,AJ438332/g' data/intermediate/blast/top_resi.csv
sed -i -e 's/EU552210/UGT2,EU552210/g' data/intermediate/blast/top_resi.csv
sed -i -e 's/AJ438330/UGT-3i11a,AJ438330/g' data/intermediate/blast/top_resi.csv
# sed -i -e 's/AJ438337/UGT-4ni13a,AJ438337/g' data/intermediate/blast/top_resi.csv # non-specific
# sed -i -e 's/AJ438334/UGT-1i15b,AJ438334/g' data/intermediate/blast/top_resi.csv # non-specific
sed -i -e 's/JX624788/UGT-CM,JX624788/g' data/intermediate/blast/top_resi.csv
sed -i -e 's/EF439840/LR1,EF439840/g' data/intermediate/blast/top_resi.csv
sed -i -e 's/AJ438333/UGT-4i13b,AJ438333/g' data/intermediate/blast/top_resi.csv
sed -i -e 's/KY825229/SR13,KY825229/g' data/intermediate/blast/top_resi.csv
sed -i -e 's/JX501669/HfrDrd,JX501669/g' data/intermediate/blast/top_resi.csv
sed -i -e 's/AF483596/Hfr1,AF483596/g' data/intermediate/blast/top_resi.csv
sed -i -e 's/FJ436983/LR34,FJ436983/g' data/intermediate/blast/top_resi.csv
sed -i -e 's/KY064064/LR22a,KY064064/g' data/intermediate/blast/top_resi.csv

# FHB1
sed -i -e 's/FN564434/FHB1-CS,FN564434/g' data/intermediate/blast/top_FHB1.csv
sed -i -e 's/KX907434.1:163305-166507/PFT,KX907434.1:163305-166507/g' data/intermediate/blast/top_FHB1.csv
sed -i -e 's/KX907434/FHB1-S3,KX907434/g' data/intermediate/blast/top_FHB1.csv

# pheno
sed -i -e 's/query,/query,accession,/g' data/intermediate/blast/top_pheno.csv
sed -i -e 's/JX877784/GLU-A3-502a,JX877784/g' data/intermediate/blast/top_pheno.csv
sed -i -e 's/AY453154/GLU-A3-a,AY453154/g' data/intermediate/blast/top_pheno.csv
sed -i -e 's/FJ549935/GLU-A3-2-1,FJ549935/g' data/intermediate/blast/top_pheno.csv
sed -i -e 's/JX877817/GLU-A3-391,JX877817/g' data/intermediate/blast/top_pheno.csv
sed -i -e 's/M22208/GLU-A1-Ax2*,M22208/g' data/intermediate/blast/top_pheno.csv
sed -i -e 's/HQ613179/GLU-A1-x,HQ613179/g' data/intermediate/blast/top_pheno.csv
sed -i -e 's/JX877811/GLU-B3-578Bb,JX877811/g' data/intermediate/blast/top_pheno.csv
sed -i -e 's/EU369722/GLU-B3-2-2,EU369722/g' data/intermediate/blast/top_pheno.csv
sed -i -e 's/JX878024/GLU-B3-593,JX878024/g' data/intermediate/blast/top_pheno.csv
sed -i -e 's/FJ755306/GLU-B3-1,FJ755306/g' data/intermediate/blast/top_pheno.csv
# sed -i -e 's/EU369731/GLU-B3-71,EU369731/g' data/intermediate/blast/top_pheno.csv # odd man out
sed -i -e 's/X13927/GLU-B1-b,X13927/g' data/intermediate/blast/top_pheno.csv
sed -i -e 's/AB263219/GLU-B1-i,AB263219/g' data/intermediate/blast/top_pheno.csv
sed -i -e 's/JF736014/GLU-B1-y8,JF736014/g' data/intermediate/blast/top_pheno.csv
sed -i -e 's/DQ357054/GLU-D3-21,DQ357054/g' data/intermediate/blast/top_pheno.csv
sed -i -e 's/DQ357057/GLU-D3-31,DQ357057/g' data/intermediate/blast/top_pheno.csv
sed -i -e 's/AB062873/GLU-D3-8-IV,AB062873/g' data/intermediate/blast/top_pheno.csv
sed -i -e 's/FJ755312/GLU-D3-5,FJ755312/g' data/intermediate/blast/top_pheno.csv
sed -i -e 's/AY585349/GLU3-S2,AY585349/g' data/intermediate/blast/top_pheno.csv
sed -i -e 's/DQ357052/GLU-D3-11,DQ357052/g' data/intermediate/blast/top_pheno.csv
sed -i -e 's/DQ457417/GLU-D3-42,DQ457417/g' data/intermediate/blast/top_pheno.csv
sed -i -e 's/EU189096/GLU-D3,EU189096/g' data/intermediate/blast/top_pheno.csv
sed -i -e 's/X12928/GLU-D1-1d,X12928/g' data/intermediate/blast/top_pheno.csv
sed -i -e 's/JX173939/GLU-D1-1,JX173939/g' data/intermediate/blast/top_pheno.csv
sed -i -e 's/BK006459/GLU-D1-12,BK006459/g' data/intermediate/blast/top_pheno.csv
sed -i -e 's/DQ885753/PPD-A1,DQ885753/g' data/intermediate/blast/top_pheno.csv
sed -i -e 's/DQ885766/PPD-D1,DQ885766/g' data/intermediate/blast/top_pheno.csv
sed -i -e 's/AY515506/PPO,AY515506/g' data/intermediate/blast/top_pheno.csv
sed -i -e 's/AB191458/Tamyb10-A1,AB191458/g' data/intermediate/blast/top_pheno.csv
sed -i -e 's/AB191459/Tamyb10-B1,AB191459/g' data/intermediate/blast/top_pheno.csv
sed -i -e 's/AB191460/Tamby10-D1,AB191460/g' data/intermediate/blast/top_pheno.csv
sed -i -e 's/JF930277/Rht-A1,JF930277/g' data/intermediate/blast/top_pheno.csv
sed -i -e 's/JX993615/Rht-B1,JX993615/g' data/intermediate/blast/top_pheno.csv
sed -i -e 's/HE585643/Rht-D1,HE585643/g' data/intermediate/blast/top_pheno.csv
sed -i -e 's/AY485974/VRN2,AY485974/g' data/intermediate/blast/top_pheno.csv
sed -i -e 's/AY747600/VRN-A1,AY747600/g' data/intermediate/blast/top_pheno.csv
sed -i -e 's/AY747604/VRN-B1,AY747604/g' data/intermediate/blast/top_pheno.csv
# sed -i -e 's/X03042/GLU1-y,X03042/g' data/intermediate/blast/top_pheno.csv # odd man out
sed -i -e 's/DQ363911/PINa,DQ363911/g' data/intermediate/blast/top_pheno.csv
sed -i -e 's/AB262660/PINb,AB262660/g' data/intermediate/blast/top_pheno.csv
sed -i -e 's/AY747606/VRN-D1,AY747606/g' data/intermediate/blast/top_pheno.csv
sed -i -e 's/DQ890162/VRN-B3,DQ890162/g' data/intermediate/blast/top_pheno.csv
sed -i -e 's/DQ885757/PPD-B1,DQ885757/g' data/intermediate/blast/top_pheno.csv