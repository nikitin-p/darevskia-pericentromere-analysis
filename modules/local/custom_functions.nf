def trimSuffix(String original, String suffix) {
    if(original.endsWith(suffix)) {
        return original.substring(0, original.length() - suffix.length())
    }
    return original
}

def extract_species(contigs_name) {
    def contigs_simple_name = contigs_name.tokenize('/')[-1]
    def m = contigs_simple_name =~ /_([NV])/
    return m[0][1]
}

def extract_reverse_species(contigs_name) {
    if (extract_species(contigs_name) == 'N') {
        return 'V'
    } else {
        return 'N'
    }
}
