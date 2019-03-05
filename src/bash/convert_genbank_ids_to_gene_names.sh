# convert ncbi identifiers to gene names
# resi
sed -i -e 's/KT962210/ACT/g' data/intermediate/blast/top_resi.csv
sed -i -e 's/KU562862/DGK/g' data/intermediate/blast/top_resi.csv
sed -i -e 's/KC244204/GLI1/g' data/intermediate/blast/top_resi.csv
sed -i -e 's/AF483596/Hfr1/g' data/intermediate/blast/top_resi.csv
sed -i -e 's/JX501669/HfrDrd/g' data/intermediate/blast/top_resi.csv
sed -i -e 's/EF439840/LR1/g' data/intermediate/blast/top_resi.csv
sed -i -e 's/U51330/LR10/g' data/intermediate/blast/top_resi.csv
sed -i -e 's/FJ876280/LR21/g' data/intermediate/blast/top_resi.csv
sed -i -e 's/KY064064/LR22a/g' data/intermediate/blast/top_resi.csv
sed -i -e 's/FJ436983/LR34/g' data/intermediate/blast/top_resi.csv
sed -i -e 's/FJ212314/PM3/g' data/intermediate/blast/top_resi.csv
sed -i -e 's/KY825229/SR13/g' data/intermediate/blast/top_resi.csv
sed -i -e 's/LN883754/SR22/g' data/intermediate/blast/top_resi.csv
sed -i -e 's/KF031303/SR33/g' data/intermediate/blast/top_resi.csv
sed -i -e 's/LN883757/SR45/g' data/intermediate/blast/top_resi.csv
sed -i -e 's/EU568801/UDP-GT/g' data/intermediate/blast/top_resi.csv
sed -i -e 's/KU562861/WRKY70/g' data/intermediate/blast/top_resi.csv
# special resi
sed -i -e 's/FN564434/FHB1/g' data/intermediate/blast/top_FHB1.csv
#pheno
sed -i -e 's/JF736016/GLU-HMW/g' data/intermediate/blast/top_pheno.csv
sed -i -e 's/AJ937920/GLU-LMW/g' data/intermediate/blast/top_pheno.csv
sed -i -e 's/DQ357055/GLU-LMW-2/g' data/intermediate/blast/top_pheno.csv
sed -i -e 's/AB181238/PINa/g' data/intermediate/blast/top_pheno.csv
sed -i -e 's/AY598029/PINb/g' data/intermediate/blast/top_pheno.csv
sed -i -e 's/DQ885753/PPD-A1/g' data/intermediate/blast/top_pheno.csv
sed -i -e 's/DQ885757/PPD-B1/g' data/intermediate/blast/top_pheno.csv
sed -i -e 's/DQ885766/PPD-D1/g' data/intermediate/blast/top_pheno.csv
sed -i -e 's/AY515506/PPO/g' data/intermediate/blast/top_pheno.csv
sed -i -e 's/JF930277/Rht-A1/g' data/intermediate/blast/top_pheno.csv
sed -i -e 's/JX993615/Rht-B1/g' data/intermediate/blast/top_pheno.csv
sed -i -e 's/HE585643/Rht-D1/g' data/intermediate/blast/top_pheno.csv
sed -i -e 's/AB191458/Tamyb10-A1/g' data/intermediate/blast/top_pheno.csv
sed -i -e 's/AB191459/Tamyb10-B1/g' data/intermediate/blast/top_pheno.csv
sed -i -e 's/AB191460/Tamby10-D1/g' data/intermediate/blast/top_pheno.csv
sed -i -e 's/AY747600/VRN-A1/g' data/intermediate/blast/top_pheno.csv
sed -i -e 's/AY747604/VRN-B1/g' data/intermediate/blast/top_pheno.csv
sed -i -e 's/DQ890162/VRN-B3/g' data/intermediate/blast/top_pheno.csv
sed -i -e 's/AY747606/VRN-D1/g' data/intermediate/blast/top_pheno.csv
sed -i -e 's/AY485974/VRN2/g' data/intermediate/blast/top_pheno.csv