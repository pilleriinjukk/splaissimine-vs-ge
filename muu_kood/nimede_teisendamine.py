tekst = "Nurlan Kerimov, James D. Hayhurst, Jonathan R. Manning, Peter Walter, Liis Kolberg, Kateryna Peikova, Marija Samoviča, Tony Burdett, Simon Jupp, Helen Parkinson, Irene Papatheodorou, Daniel R. Zerbino, Kaur Alasoo"
nimed = tekst.split(", ")
tulemus = ""
for nimi in nimed:
    osad = nimi.split(" ")
    print(osad[-1])
    tulemus += osad[-1] + ", "
    for osa in osad[:-1]:
        tulemus += osa[0] + ". "
    tulemus = tulemus[:-1]
    tulemus += ", "
print(tulemus)
