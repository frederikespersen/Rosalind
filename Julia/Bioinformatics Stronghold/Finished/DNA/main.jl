function count_nbs(dna::String)
    nbs = Dict{Char, Int}(
        'A' => 0,
        'C' => 0,
        'G' => 0,
        'T' => 0)
    
    for nb in dna
        nbs[nb] += 1
    end

    return nbs

end

function load_file(filepath::String)
    cd("Julia/Bioinformatics Stronghold/DNA")
    file = open(filepath)
    lines = readlines(file)
    return lines
end


dna = load_file("sequence.txt")[1]
counts = count_nbs(dna)
counts = [counts['A'], counts['C'], counts['G'], counts['T']]
counts = map(x->string(x), counts)
print(join(counts, ' '))