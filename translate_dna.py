from re import A

dna = input('enter a DNA template: ')
#'TTAACTGACGCTTGACATAGCTGATG'
 #caucagcuAUGUCAAGCGUCAGUUAA

def comp_strand(query):
    query = query.upper()
    comp = ''
    for letter in query:
        if letter == 'A':
            comp += 'T'
        elif letter == 'T':
            comp += 'A'
        elif letter == 'C':
            comp += 'G'
        elif letter == 'G':
            comp += 'C'
    print("       3' end " + comp + " 5' end")
    return (comp[::-1])

def transcribe(template):
    coding = template.replace('T', 'U')
    rna = ''; b=0
    exon = False
    while(b < len(coding)):
        if coding[b:b+3] == 'AUG':
            exon = True
        if coding[b:b+3] == 'UAG' or coding[b:b+3] == 'UGA' or coding[b:b+3] == 'UAA':
            exon = False
        if exon == True:
            rna += coding[b:b+3].upper()
            b+=3
        if exon == False:
            rna += coding[b].lower()
            b+=1
    mrna = ''
    for i in range(len(rna)-2):
        if (rna[i:i+3] == 'uag' or rna[i:i+3] == 'uga' or rna[i:i+3] == 'uaa') and rna[i-1].isupper():
            mrna += rna[i:i+3].upper()
            i+=2
        else:
            mrna += rna[i]
    return mrna

def splice(rna):
    mrna = ''
    for bp in rna:
        if bp.isupper():
            mrna += bp
    return mrna

def translate(rna):
    codons = ''

    for nt in range(0,len(rna)-3,3):
        if rna[nt] == 'A':
            if rna[nt+1:nt+3] == 'AA' or rna[nt+1:nt+3] == 'AG': codons += ('Lys')
            elif rna[nt+1:nt+3] == 'AC' or rna[nt+1:nt+3] == 'AU': codons += ('Asn')
            elif rna[nt+1:nt+3] == 'GA' or rna[nt+1:nt+3] == 'GG': codons += ('Arg');
            elif rna[nt+1:nt+3] == 'GC' or rna[nt+1:nt+3] == 'GU': codons += ('Ser')
            elif rna[nt+1] == 'C': codons += ('Thr')
            elif rna[nt+1:nt+3] == 'UG': codons += ('MET');
            elif rna[nt+1:nt+3] == 'UA': codons += ('Ile')
            elif rna[nt+1:nt+3] == 'UC': codons += ('Ile')
            elif rna[nt+1:nt+3] == 'UU': codons += ('Ile')
        elif rna[nt] == 'G':
            if rna[nt+1:nt+3] == 'AA' or rna[nt+1:nt+3] == 'AG': codons += ('Glu')
            elif rna[nt+1:nt+3] == 'AC' or rna[nt+1:nt+3] == 'AU': codons += ('Asp')
            elif rna[nt+1] == 'G': codons += ('Gly')
            elif rna[nt+1] == 'C': codons += ('Ala')
            elif rna[nt+1] == 'U': codons += ('Val')
        elif rna[nt] == 'C':
            if rna[nt+1:nt+3] == 'AA' or rna[nt+1:nt+3] == 'AG': codons += ('Gln')
            elif rna[nt+1:nt+3] == 'AC' or rna[nt+1:nt+3] == 'AU': codons += ('His')
            elif rna[nt+1] == 'G': codons += ('Arg')
            elif rna[nt+1] == 'C': codons += ('Pro')
            elif rna[nt+1] == 'U': codons += ('Leu')
        elif rna[nt] == 'U':
            if rna[nt+1] == 'C': codons += ('Ser')
            elif rna[nt+1:nt+3] == 'AA' or rna[nt+1:nt+3] == 'AG':
                codons += ('STOP\n');
            elif rna[nt+1:nt+3] == 'AC' or rna[nt+1:nt+3] == 'AU': codons += ('Tyr')
            if rna[nt+1:nt+3] == 'GA':
                codons += ('STOP\n');
            elif rna[nt+1:nt+3] == 'GG': codons += ('Trp')
            elif rna[nt+1:nt+3] == 'GC' or rna[nt+1:nt+3] == 'GU': codons += ('Cys')
            elif rna[nt+1:nt+3] == 'UA' or rna[nt+1:nt+3] == 'UG': codons += ('Leu')
            elif rna[nt+1:nt+3] == 'UC' or rna[nt+1:nt+3] == 'UU': codons += ('Phe')
    return codons

print("  Template 5' "+ dna + " 3' end")

complement = comp_strand(dna)
print("    Coding 5' " + complement + " 3' end")

mRNA = transcribe(complement)
print('Transcipt seq ' + mRNA)

mature = splice(mRNA)
print('   Coding seq ' + mature)

proteins = translate(mature)
print('  Protein seq ' + proteins)