function load_string(rel_filepath::String)
    filepath = string(pwd(), '/', rel_filepath)
    file = open(filepath)
    s = read(file, String)
    return s
end

function complimentary(dna::String)
    comps = Dict{Char, Char}(
        'A'=>'T',
        'T'=>'A',
        'C'=>'G',
        'G'=>'C')
    compliment_rev = map(x->comps[x], dna)
    compliment = reverse(compliment_rev)
    return compliment
end

dna = load_string("Julia/Bioinformatics Stronghold/REVC/dna.txt")
comp = complimentary(dna)

println(comp)
