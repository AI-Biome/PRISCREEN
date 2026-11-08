from taxonomy import lowest_common_rank, parse_header_taxonomy, rank_index


def test_panel_header():
    t = parse_header_taxonomy("loc1|species=Yersinia_pestis|genome=Yersinia_pestis_2")
    assert t["species"] == "Yersinia pestis" and t["genus"] == "Yersinia" and t["domain"] is None


def test_panel_single_token_species_has_no_genus():
    t = parse_header_taxonomy("amp1|species=SpeciesA")
    assert t["species"] == "SpeciesA" and t["genus"] is None


def test_silva_seven_levels():
    t = parse_header_taxonomy("AB001.1.1500 Bacteria;Pseudomonadota;Gammaproteobacteria;Enterobacterales;"
                              "Enterobacteriaceae;Escherichia-Shigella;Escherichia coli")
    assert t["domain"] == "Bacteria" and t["family"] == "Enterobacteriaceae"
    assert t["genus"] == "Escherichia-Shigella" and t["species"] == "Escherichia coli"


def test_silva_long_lineage_uses_positions_from_the_end():
    t = parse_header_taxonomy("AB2.1.1400 Eukaryota;SAR;Alveolata;Dinoflagellata;Dinophyceae;Gonyaulacales;"
                              "Gonyaulacaceae;Alexandrium;Alexandrium minutum")
    assert t["domain"] == "Eukaryota" and t["genus"] == "Alexandrium" and t["species"] == "Alexandrium minutum"
    assert t["phylum"] is None


def test_gtdb_prefixed_with_annotations():
    t = parse_header_taxonomy("RS_GCF_000006945.2~NC_003197.2 d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;"
                              "o__Enterobacterales;f__Enterobacteriaceae;g__Salmonella;s__Salmonella enterica [location=1..1500]")
    assert t["genus"] == "Salmonella" and t["species"] == "Salmonella enterica" and t["phylum"] == "Pseudomonadota"


def test_sintax_style():
    t = parse_header_taxonomy("X1;tax=d:Bacteria,p:Firmicutes,c:Bacilli,o:Bacillales,f:Bacillaceae,g:Bacillus,s:Bacillus_subtilis;")
    assert t["species"] == "Bacillus subtilis" and t["order"] == "Bacillales"


def test_uncultured_is_unknown():
    t = parse_header_taxonomy("A Bacteria;Bacillota;Bacilli;Bacillales;Bacillaceae;Bacillus;uncultured bacterium")
    assert t["species"] is None and t["genus"] == "Bacillus"


def test_plain_id_becomes_species():
    assert parse_header_taxonomy("plainid")["species"] == "plainid"


def test_lca_skips_ranks_unknown_everywhere():
    a, b = parse_header_taxonomy("l|species=Yersinia_pestis"), parse_header_taxonomy("l|species=Yersinia_pseudotuberculosis")
    assert lowest_common_rank([a]) == ("species", "Yersinia pestis")
    assert lowest_common_rank([a, b]) == ("genus", "Yersinia")
    assert lowest_common_rank([parse_header_taxonomy("a|species=SpA"), parse_header_taxonomy("a|species=SpB")]) == (None, None)


def test_lca_stops_when_some_are_unknown():
    a = parse_header_taxonomy("A Bacteria;P;C;O;F;Bacillus;Bacillus subtilis")
    b = parse_header_taxonomy("B Bacteria;P;C;O;F;Bacillus;uncultured bacterium")
    assert lowest_common_rank([a, b]) == ("genus", "Bacillus")


def test_rank_index():
    assert rank_index("genus") == 5 and rank_index("kingdom") == 0 and rank_index("nonsense") == -1


def test_silva_short_lineage_keeps_organism_as_species():
    t = parse_header_taxonomy("X2 Bacteria;Bacteroidota;Bacteroidia;Bacteroidales;Alistipes finegoldii")
    assert t["species"] == "Alistipes finegoldii" and t["genus"] == "Alistipes" and t["family"] is None
    assert t["order"] == "Bacteroidales" and t["domain"] == "Bacteria"
    assert parse_header_taxonomy("X3 Bacteria;Alistipes sp.")["species"] == "Alistipes sp."


def test_pipe_in_id_with_description_lineage():
    t = parse_header_taxonomy("RS_GCF_1|NC_1 d__Bacteria;p__P;c__C;o__O;f__F;g__Salmonella;s__Salmonella enterica [location=1..9]")
    assert t["species"] == "Salmonella enterica" and t["domain"] == "Bacteria"
    t = parse_header_taxonomy("amp1|species=Yersinia_pestis|genome=Y_1 accn|CP1:1-2 strand=+ size=2")
    assert t["species"] == "Yersinia pestis" and t["genus"] == "Yersinia"
