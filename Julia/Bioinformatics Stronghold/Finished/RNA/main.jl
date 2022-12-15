function load_string(rel_filepath::String)
    filepath = string(pwd(), '/', rel_filepath)
    file = open(filepath)
    s = read(file, String)
    return s
end

function translate(dna)
    rna = replace(dna, 'T'=>'U')
    return rna
end

dna = load_string("Julia/Bioinformatics Stronghold/RNA/dna.txt")
rna = translate(dna)

print(rna)