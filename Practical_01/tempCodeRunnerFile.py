from Bio import SeqIO
import csv
import matplotlib.pyplot as plt

# Enter the path to your downloaded GenBank file
file_path = "C:\Gitfolder\Practical_01\genomic.gbff"

# Genes related to plant growth promotion
target_genes = [
    "nifH", "nifD", "nifK", "nifA",
    "gcd", "pqqC", "pqqE", "ipdC"
]

results = []

# Read annotated genome file
for record in SeqIO.parse(file_path, "genbank"):
    for feature in record.features:

        if feature.type not in ["gene", "CDS"]:
            continue

        qualifiers = feature.qualifiers

        gene = qualifiers.get("gene", [""])[0]
        product = qualifiers.get("product", ["Unknown"])[0]
        locus = qualifiers.get("locus_tag", ["Unknown"])[0]

        for target in target_genes:
            if gene.lower() == target.lower():
                results.append({
                    "Gene": target,
                    "Product": product,
                    "Locus_tag": locus,
                    "Start": int(feature.location.start) + 1,
                    "End": int(feature.location.end)
                })

# Save results to CSV
with open("gene_results.csv", "w", newline="",
          encoding="utf-8") as file:

    columns = [
        "Gene", "Product", "Locus_tag", "Start", "End"
    ]

    writer = csv.DictWriter(file, fieldnames=columns)
    writer.writeheader()
    writer.writerows(results)

print("Analysis completed!")
print("Matching annotated features:", len(results))
print("Results saved in gene_results.csv")

# Plot gene counts
counts = {
    gene: sum(1 for row in results if row["Gene"] == gene)
    for gene in target_genes
}

plt.figure(figsize=(10, 5))
plt.bar(counts.keys(), counts.values())
plt.xlabel("Target gene")
plt.ylabel("Number of matching annotated features")
plt.title("PlantGrowth Gene Analyzer")
plt.tight_layout()
plt.savefig("gene_comparison.png")
plt.show()
