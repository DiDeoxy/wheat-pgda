# convert ncbi identifiers to gene names
# resi
sed -i -e 's/AY325736/PM3-B/g' data/intermediate/blast/top_resi.csv
sed -i -e 's/KF572030/PM8/g' data/intermediate/blast/top_resi.csv
sed -i -e 's/MH077963/PM17/g' data/intermediate/blast/top_resi.csv
sed -i -e 's/FJ876280/LR21/g' data/intermediate/blast/top_resi.csv
sed -i -e 's/AY270159/RGA2/g' data/intermediate/blast/top_resi.csv #~LR10
sed -i -e 's/KF031303/SR33/g' data/intermediate/blast/top_resi.csv
sed -i -e 's/LN883757/SR45/g' data/intermediate/blast/top_resi.csv
sed -i -e 's/HM133634/UGT-S3/g' data/intermediate/blast/top_resi.csv
sed -i -e 's/AK455349/UGT-CS/g' data/intermediate/blast/top_resi.csv
sed -i -e 's/KU562861/WRKY70/g' data/intermediate/blast/top_resi.csv
sed -i -e 's/KT962210/ACT/g' data/intermediate/blast/top_resi.csv
sed -i -e 's/KC244204/GLI1/g' data/intermediate/blast/top_resi.csv
sed -i -e 's/JN631792/YR5/g' data/intermediate/blast/top_resi.csv
sed -i -e 's/AJ438338/UGT-1i15e/g' data/intermediate/blast/top_resi.csv
sed -i -e 's/AJ438331/UGT-3i11b/g' data/intermediate/blast/top_resi.csv
sed -i -e 's/FJ535237/UGT-WS-2/g' data/intermediate/blast/top_resi.csv
sed -i -e 's/KX274091/UGT-UNKNOWN/g' data/intermediate/blast/top_resi.csv
sed -i -e 's/EU568801/UGT1/g' data/intermediate/blast/top_resi.csv
sed -i -e 's/AJ438332/UGT-3i11e/g' data/intermediate/blast/top_resi.csv
sed -i -e 's/EU552210/UGT2/g' data/intermediate/blast/top_resi.csv
sed -i -e 's/AJ438330/UGT-3i11a/g' data/intermediate/blast/top_resi.csv
# sed -i -e 's/AJ438337/UGT-4ni13a/g' data/intermediate/blast/top_resi.csv # non-specific
# sed -i -e 's/AJ438334/UGT-1i15b/g' data/intermediate/blast/top_resi.csv # non-specific
sed -i -e 's/JX624788/UGT-CM/g' data/intermediate/blast/top_resi.csv
sed -i -e 's/EF439840/LR1/g' data/intermediate/blast/top_resi.csv
sed -i -e 's/AJ438333/UGT-4i13b/g' data/intermediate/blast/top_resi.csv
sed -i -e 's/KY825229/SR13/g' data/intermediate/blast/top_resi.csv
sed -i -e 's/JX501669/HfrDrd/g' data/intermediate/blast/top_resi.csv
sed -i -e 's/AF483596/Hfr1/g' data/intermediate/blast/top_resi.csv
sed -i -e 's/FJ436983/LR34/g' data/intermediate/blast/top_resi.csv
sed -i -e 's/KY064064/LR22a/g' data/intermediate/blast/top_resi.csv


# FHB1
sed -i -e 's/FN564434/FHB1-CS/g' data/intermediate/blast/top_FHB1.csv
sed -i -e 's/KX907434.1:163305-166507/PFT/g' data/intermediate/blast/top_FHB1.csv
sed -i -e 's/KX907434/FHB1-S3/g' data/intermediate/blast/top_FHB1.csv

# pheno
sed -i -e 's/JX877784/GLU-A3-502a/g' data/intermediate/blast/top_pheno.csv
sed -i -e 's/AY453154/GLU-A3-a/g' data/intermediate/blast/top_pheno.csv
sed -i -e 's/FJ549935/GLU-A3-2-1/g' data/intermediate/blast/top_pheno.csv
sed -i -e 's/JX877817/GLU-A3-391/g' data/intermediate/blast/top_pheno.csv
sed -i -e 's/M22208/GLU-A1-Ax2*/g' data/intermediate/blast/top_pheno.csv
sed -i -e 's/HQ613179/GLU-A1-x/g' data/intermediate/blast/top_pheno.csv
sed -i -e 's/JX877811/GLU-B3-578Bb/g' data/intermediate/blast/top_pheno.csv
sed -i -e 's/EU369722/GLU-B3-2-2/g' data/intermediate/blast/top_pheno.csv
sed -i -e 's/JX878024/GLU-B3-593/g' data/intermediate/blast/top_pheno.csv
sed -i -e 's/FJ755306/GLU-B3-1/g' data/intermediate/blast/top_pheno.csv
# sed -i -e 's/EU369731/GLU-B3-71/g' data/intermediate/blast/top_pheno.csv # odd man out
sed -i -e 's/X13927/GLU-B1-b/g' data/intermediate/blast/top_pheno.csv
sed -i -e 's/AB263219/GLU-B1-i/g' data/intermediate/blast/top_pheno.csv
sed -i -e 's/JF736014/GLU-B1-y8/g' data/intermediate/blast/top_pheno.csv
sed -i -e 's/DQ357054/GLU-D3-21/g' data/intermediate/blast/top_pheno.csv
sed -i -e 's/DQ357057/GLU-D3-31/g' data/intermediate/blast/top_pheno.csv
sed -i -e 's/AB062873/GLU-D3-8-IV/g' data/intermediate/blast/top_pheno.csv
sed -i -e 's/FJ755312/GLU-D3-5/g' data/intermediate/blast/top_pheno.csv
sed -i -e 's/AY585349/GLU3-S2/g' data/intermediate/blast/top_pheno.csv
sed -i -e 's/DQ357052/GLU-D3-11/g' data/intermediate/blast/top_pheno.csv
sed -i -e 's/DQ457417/GLU-D3-42/g' data/intermediate/blast/top_pheno.csv
sed -i -e 's/EU189096/GLU-D3/g' data/intermediate/blast/top_pheno.csv
sed -i -e 's/X12928/GLU-D1-1d/g' data/intermediate/blast/top_pheno.csv
sed -i -e 's/JX173939/GLU-D1-1/g' data/intermediate/blast/top_pheno.csv
sed -i -e 's/BK006459/GLU-D1-12/g' data/intermediate/blast/top_pheno.csv
sed -i -e 's/DQ885753/PPD-A1/g' data/intermediate/blast/top_pheno.csv
sed -i -e 's/DQ885766/PPD-D1/g' data/intermediate/blast/top_pheno.csv
sed -i -e 's/AY515506/PPO/g' data/intermediate/blast/top_pheno.csv
sed -i -e 's/AB191458/Tamyb10-A1/g' data/intermediate/blast/top_pheno.csv
sed -i -e 's/AB191459/Tamyb10-B1/g' data/intermediate/blast/top_pheno.csv
sed -i -e 's/AB191460/Tamby10-D1/g' data/intermediate/blast/top_pheno.csv
sed -i -e 's/JF930277/Rht-A1/g' data/intermediate/blast/top_pheno.csv
sed -i -e 's/JX993615/Rht-B1/g' data/intermediate/blast/top_pheno.csv
sed -i -e 's/HE585643/Rht-D1/g' data/intermediate/blast/top_pheno.csv
sed -i -e 's/AY485974/VRN2/g' data/intermediate/blast/top_pheno.csv
sed -i -e 's/AY747600/VRN-A1/g' data/intermediate/blast/top_pheno.csv
sed -i -e 's/AY747604/VRN-B1/g' data/intermediate/blast/top_pheno.csv
# sed -i -e 's/X03042/GLU1-y/g' data/intermediate/blast/top_pheno.csv # odd man out
sed -i -e 's/AB181238/PINa/g' data/intermediate/blast/top_pheno.csv
sed -i -e 's/AY598029/PINb/g' data/intermediate/blast/top_pheno.csv
sed -i -e 's/AY747606/VRN-D1/g' data/intermediate/blast/top_pheno.csv
sed -i -e 's/DQ890162/VRN-B3/g' data/intermediate/blast/top_pheno.csv
sed -i -e 's/DQ885757/PPD-B1/g' data/intermediate/blast/top_pheno.csv