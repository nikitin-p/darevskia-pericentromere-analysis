def trimSuffix(String original, String suffix) {
    if(original.endsWith(suffix)) {
        return original.substring(0, original.length() - suffix.length())
    }
    return original
}

def extract_species(contigs_name) {
    def m = contigs_name =~ /contigs_([NV])/;
    return m[0][1]​
}