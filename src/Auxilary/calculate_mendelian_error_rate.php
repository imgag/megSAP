<?php

error_reporting(E_ERROR | E_WARNING | E_PARSE | E_NOTICE);

require_once(dirname($_SERVER['SCRIPT_FILENAME'])."/../Common/all.php");

// parse command line arguments
$parser = new ToolBase("calculate_mendelian_error_rate", "Calculates the mendelian error rate of a given trio VCF file.");
$parser->addInfile("in",  "Input trio VCF file.", false);
extract($parser->parse($argv));


// Function to extract genotype (GT field) from a FORMAT column
function extract_genotype($sample_field)
{
    $fields = explode(":", $sample_field);
    return isset($fields[0]) ? $fields[0] : null;
}

// Function to check if the child's genotype follows Mendelian inheritance
function is_mendelian($child, $parent1, $parent2)
{
    // Split genotypes into alleles
    $child_alleles = explode("/", $child);
    $parent1_alleles = explode("/", $parent1);
    $parent2_alleles = explode("/", $parent2);

    // Mendelian inheritance: child's alleles must be a combination of one allele from each parent
    foreach ($child_alleles as $allele) 
    {
        if (!in_array($allele, $parent1_alleles) && !in_array($allele, $parent2_alleles)) 
        {
            return false; // Error found
        }
    }
    return true;
}

// Initialize counters for Mendelian errors and total variants
$total_variants = 0;
$mendelian_errors = 0;

// Open the VCF file for reading
if (($handle = fopen($in, "r")) !== false) 
{
    while (($line = fgets($handle)) !== false)
    {
        // Skip header lines
        if (substr($line, 0, 1) === "#")
        {
            continue;
        }

        // Split the line into columns
        $columns = explode("\t", trim($line));

        // Extract relevant columns //TODO add child sample name as parameter and get column ID
        if (empty($columns[9])) echo $line; 
        $child_genotype = extract_genotype($columns[11]);  // NA12878x2_80
        $parent1_genotype = extract_genotype($columns[10]);  // NA12891_14
        $parent2_genotype = extract_genotype($columns[9]);  // NA12892_18

        // If genotypes are missing, skip the variant
        if ($child_genotype === null || $parent1_genotype === null || $parent2_genotype === null)
        {
            continue;
        }

        $total_variants++;

        // Check for Mendelian error
        if (!is_mendelian($child_genotype, $parent1_genotype, $parent2_genotype)) 
        {
            $mendelian_errors++;
        }
    }
    fclose($handle);
}

// Calculate and print the Mendelian error rate
if ($total_variants > 0) 
{
    $error_rate = $mendelian_errors / $total_variants;
    echo "Mendelian Error Rate: " . ($error_rate * 100) . "%\n";
    echo "Total Variants: $total_variants\n";
    echo "Mendelian Errors: $mendelian_errors\n";
}
else 
{
    echo "No valid variants found.\n";
}

?>