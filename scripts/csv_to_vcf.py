import csv
import argparse

def csv_to_vcf(input_csv, output_vcf):
    try:
        # D'abord détecter le séparateur
        with open(input_csv, 'r') as f:
            dialect = csv.Sniffer().sniff(f.read(1024))
            f.seek(0)
            
            # Lire la première ligne pour vérifier les colonnes
            reader = csv.DictReader(f, delimiter=dialect.delimiter)
            fieldnames = reader.fieldnames
            
            # Vérifier les colonnes nécessaires (avec différentes orthographes possibles)
            required_columns = {'pos', 'refcall', 'alt', 'af', 'depth', 'heteroplasmy', 'type', 'allele_count'}
            available_columns = {col.lower() for col in fieldnames}
            
            if not required_columns.issubset(available_columns):
                print("Erreur: Colonnes manquantes dans le fichier CSV")
                print(f"Colonnes requises: {required_columns}")
                print(f"Colonnes disponibles: {available_columns}")
                return False
            
            # Réouvrir le fichier pour traitement complet
            f.seek(0)
            reader = csv.DictReader(f, delimiter=dialect.delimiter)
            
            # Écrire l'en-tête VCF
            with open(output_vcf, 'w') as vcffile:
                vcf_header = """##fileformat=VCFv4.2
##source=csv2vcf_converter
##reference=rCRS
##INFO=<ID=AF,Number=1,Type=Float,Description="Allele Frequency">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##INFO=<ID=HT,Number=1,Type=Float,Description="Heteroplasmy">
##INFO=<ID=AC,Number=1,Type=Integer,Description="Allele Count">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE
"""
                vcffile.write(vcf_header)
                
                for row in reader:
                    # Trouver les noms de colonnes réels (insensibles à la casse)
                    col_map = {col.lower(): col for col in fieldnames}
                    
                    chrom = 'MT'
                    pos = row[col_map['pos']]
                    ref = row[col_map['refcall']]
                    alt = row[col_map['alt']]
                    af = row[col_map['af']]
                    depth = row[col_map['depth']]
                    heteroplasmy = row[col_map['heteroplasmy']]
                    ac = row[col_map['allele_count']]
                    var_type = row[col_map['type']].lower()
                    
                    info = f"AF={af};DP={depth};HT={heteroplasmy};AC={ac}"
                    
                    if var_type == 'insertion':
                        alt = ref + alt
                    elif var_type == 'deletion':
                        alt = ref
                        ref = ref + row[col_map['alt']]
                    
                    vcf_line = f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t.\t{info}\tGT\t0/1\n"
                    vcffile.write(vcf_line)
            
            return True
            
    except Exception as e:
        print(f"Erreur lors de la conversion: {str(e)}")
        return False

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Convertir un fichier CSV de variants mitochondriaux en format VCF pour Haplogrep')
    parser.add_argument('-i', '--input', required=True, help='Fichier CSV d\'entrée')
    parser.add_argument('-o', '--output', required=True, help='Fichier VCF de sortie')
    
    args = parser.parse_args()
    
    print(f"Conversion de {args.input} en {args.output}...")
    if csv_to_vcf(args.input, args.output):
        print("Conversion terminée avec succès!")
    else:
        print("La conversion a échoué. Voir les messages d'erreur ci-dessus.")
