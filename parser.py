def load_annotations(data_folder):
    import sys
    import re
    import xml.etree.ElementTree as ET

    open(data_folder + "/rhea.rdf", "r")
    tree = ET.parse(rhea_rdf)
    root = tree.getroot()

    #dictionary with namespaces
    ns = {"rh":"http://rdf.rhea-db.org/",
          "rdfs":"http://www.w3.org/2000/01/rdf-schema#",
          "rdf":"http://www.w3.org/1999/02/22-rdf-syntax-ns#"}
    
    #Loop for assembling compound_lib
    def rp_name_func(description):
        node = description.find("rh:name", ns)
        if node is None:
            pass
        else:
            rp_name = node.text
            return rp_name


    def rp_formula_func(description):
        node = description.find("rh:formula", ns)
        if node is None:
            pass
        else:
            rp_formula = node.text
            return rp_formula


    def rp_charge_func(description):
        node = description.find("rh:charge", ns)
        if node is None:
            pass
        else:
            rp_charge = node.text
            return rp_charge


    def rp_chebi_func(description):
        node = description.find("rh:chebi", ns)
        if node is None:
            pass
        else:
            rp_chebi = node.attrib["{http://www.w3.org/1999/02/22-rdf-syntax-ns#}resource"].lstrip("http://purl.obolibrary.org/obo/").replace("_", ":")
            return rp_chebi


    def add_reactive_part(description):
        comp_num = re.sub(pattern = r"_rp\d", repl = "", string = rp).lstrip("Compound_")
        comp_index = None
        for compound in compound_lib:
            if comp_num == compound.get("comp_num"):
                comp_index = compound_lib.index(compound)
        if comp_index is None:
            print("Assigning reactive part to compound entry not made yet")
        else:
            reactive_part = {}
            rp_name = rp_name_func(description)
            rp_formula = rp_formula_func(description)
            rp_charge = rp_charge_func(description)
            rp_chebi = rp_chebi_func(description)
            reactive_part["name"] = rp_name
            reactive_part["formula"] = rp_formula
            reactive_part["charge"] = rp_charge
            reactive_part["chebi_id"] = rp_chebi
            compound_lib[comp_index]["reactive_parts"].append(reactive_part)


    def compound_name(comp_entry):
        node = description.find("rh:name", ns)
        if node is None:
            pass
        else:
            compound_name = node.text
            comp_entry["name"] = compound_name


    def compound_formula(comp_entry):
        node = description.find("rh:formula", ns)
        if node is None:
            pass
        else:
            formula = node.text
            if formula is None:
                pass
            else:
                formula = formula.rstrip("<i><sub>n</sub></i>")
                comp_entry["formula"] = formula


    def compound_charge(comp_entry):
        node = description.find("rh:charge", ns)
        if node is None:
            pass
        else:
            compound_charge = node.text.rstrip("<i><sub>n</sub></i>")
            comp_entry["charge"] = compound_charge


    #Fills chebi library with chebi entries
    def gather_chebi_data(comp_num):
        comp_entry = {}
        comp_entry["comp_num"] = comp_num
        chebi_id = description.find("rh:accession", ns).text
        comp_entry["chebi_id"] = chebi_id
        compound_name(comp_entry)
        compound_formula(comp_entry)
        compound_charge(comp_entry)
        return comp_entry        


    def gather_generic_data(comp_num):
        generic_entry = {}
        generic_entry["comp_num"] = comp_num
        generic_id = description.find("rh:accession", ns).text.lstrip("GENERIC:")
        generic_entry["generic_id"] = generic_id
        compound_name(generic_entry)
        compound_formula(generic_entry)
        compound_charge(generic_entry)
        generic_entry["reactive_parts"] = []
        return generic_entry


    def gather_poly_data(comp_num):
        poly_entry = {}
        poly_entry["comp_num"] = comp_num
        poly_id = description.find("rh:accession", ns).text.lstrip("POLYMER:")
        poly_entry["poly_id"] = poly_id
        compound_name(poly_entry)
        compound_formula(poly_entry)
        compound_charge(poly_entry)
        return poly_entry    


    compound_lib = []
    for description in root.findall("rdf:Description", ns):
        accession_node = description.find("rh:accession", ns)
        if accession_node is None:
            rp = description.attrib["{http://www.w3.org/1999/02/22-rdf-syntax-ns#}about"].lstrip("http://rdf.rhea-db.org/")
            pattern2 = re.compile("Compound_\d+_rp\d")
            if pattern2.match(rp):
                add_reactive_part(description)
            else:
                continue
        else:
            accession = accession_node.text
            if "CHEBI:" in accession:
                comp_num = description.attrib["{http://www.w3.org/1999/02/22-rdf-syntax-ns#}about"].lstrip("http://rdf.rhea-db.org/Compound_")
                comp_entry = gather_chebi_data(comp_num)
                compound_lib.append(comp_entry)

            elif "GENERIC:" in accession:
                comp_num = description.attrib["{http://www.w3.org/1999/02/22-rdf-syntax-ns#}about"].lstrip("http://rdf.rhea-db.org/Compound_")
                generic_entry = gather_generic_data(comp_num)
                compound_lib.append(generic_entry)

            elif "POLYMER:" in accession:
                comp_num = description.attrib["{http://www.w3.org/1999/02/22-rdf-syntax-ns#}about"].lstrip("http://rdf.rhea-db.org/Compound_")
                poly_entry = gather_poly_data(comp_num)
                compound_lib.append(poly_entry)

            else:
                continue   
            
    #Loop for assembling rhea_lib

    #Produces a list of all children tags of the one description
    def gather_children_tags(description):
        children_tags = []
        children = description.findall("*")
        for child in children:
            children_tags.append(child.tag)
        return children_tags


    def if_transport(comp_num, compound_entry):
        location_in = "in" in comp_num
        location_out = "out" in comp_num
        if location_in:
            compound_entry["location"] = "in"
        elif location_out:
            compound_entry["location"] = "out"
        else:
            pass


    def add_side(rhea_index, compound_num, stoich, side):
        comp_num = compound_num.group(1)
        rhea_lib[rhea_index][side].append({})
        for compound_entry in reversed(rhea_lib[rhea_index][side]):
            if len(compound_entry) == 0:
                compound_entry["stoich"] = stoich
                if_transport(comp_num, compound_entry)
                for comp in reversed(compound_lib):
                    if comp_num.rstrip("_in").rstrip("_out").lstrip("compound_") == comp.get("comp_num"):
                        compound_entry.update(comp)


    #used with descriptions of reaction sides to extract compound numbers and stoichiometry
    def gather_participant_data(description):
        children_tags = gather_children_tags(description)
        pattern = "{http://rdf.rhea-db.org/}contains[1-9]"
        for child in children_tags:
            if re.match(pattern, child):
                #Grabbing compound information
                compound_value = description.find(child).attrib["{http://www.w3.org/1999/02/22-rdf-syntax-ns#}resource"]
                pattern = ".+(compound_.*)"
                compound_num = re.match(pattern, compound_value)
                stoich = child[-1]

                #Adding to the proper rhea entry
                rhea_id = rhea_reaction_side.lstrip("http://rdf.rhea-db.org/").rstrip("_R").rstrip("_L")
                rhea_id = "RHEA:" + rhea_id
                for rhea in reversed(rhea_lib):
                    if rhea_id == rhea.get("rhea_id"):
                        rhea_index = rhea_lib.index(rhea)

                if rmc == "L":
                    add_side(rhea_index, compound_num, stoich, side="side_l")
                if rmc == "R":
                    add_side(rhea_index, compound_num, stoich, side="side_r")
                else:
                    continue


    def rhea_equation(rhea_entry):
        node = description.find("rh:equation", ns)
        if node is None:
            pass
        else:
            rhea_entry["equation"] = node.text


    def is_transport(rhea_entry):
        node = description.find("rh:isTransport", ns)
        if node is None:
            pass
        else:
            rhea_entry["is_transport"] = node.text


    def ec_link(rhea_entry):
        node = description.find("rh:ec", ns)
        if node is None:
            pass
        else:
            rhea_entry["ec_link"] = node.attrib["{http://www.w3.org/1999/02/22-rdf-syntax-ns#}resource"]
            rhea_entry["ec_id"] = rhea_entry["ec_link"].lstrip("http://purl.uniprot.org/enzyme/")    


    def status(rhea_entry):
        node = description.find("rh:status", ns)
        if node is None:
            pass
        else:
            rhea_entry["status"] = node.attrib["{http://www.w3.org/1999/02/22-rdf-syntax-ns#}resource"].lstrip("http://rdf.rhea-db.org/")


    def citations(rhea_entry):
        node = description.findall("rh:citation", ns)
        if bool(node) is False:
            pass
        else:
            rhea_entry["citations"] = []
            citations = node
            for citation in citations:
                citation = citation.attrib["{http://www.w3.org/1999/02/22-rdf-syntax-ns#}resource"].lstrip("http://rdf.ncbi.nlm.nih.gov/pubmed/")
                citation = "PMID:" + citation
                rhea_entry["citations"].append(citation)                   


    def gather_children_rheas(rhea_entry):
        if description.find("rh:directionalReaction", ns) is None:
            pass
        else:
            rhea_entry["children_rheas"] = []
            child_rheas = description.findall("rh:directionalReaction", ns)
            for child_rhea in child_rheas:
                child_rhea_name = "RHEA:" + child_rhea.attrib["{http://www.w3.org/1999/02/22-rdf-syntax-ns#}resource"].lstrip("http://rdf.rhea-db.org/")
                rhea_entry["children_rheas"].append(child_rhea_name)
            child_rhea = description.find("rh:bidirectionalReaction", ns)
            child_rhea_name = "RHEA:" + child_rhea.attrib["{http://www.w3.org/1999/02/22-rdf-syntax-ns#}resource"].lstrip("http://rdf.rhea-db.org/")
            rhea_entry["children_rheas"].append(child_rhea_name)        


    #Fills rhea dictionaries
    def gather_rhea_data(rhea_id):
        rhea_entry = {}
        rhea_entry["rhea_id"] = rhea_id
        rhea_equation(rhea_entry)
        is_transport(rhea_entry)
        ec_link(rhea_entry)
        status(rhea_entry)
        citations(rhea_entry)
        rhea_entry["side_l"] = []
        rhea_entry["side_r"] = []
        gather_children_rheas(rhea_entry)
        return rhea_entry                


    rhea_lib = []   
    for description in root.findall("rdf:Description", ns):
        node = description.find("rh:accession", ns)
        if node is None:
            rhea_reaction_side = description.attrib["{http://www.w3.org/1999/02/22-rdf-syntax-ns#}about"]
            #Right most character denoting reaction side
            rmc = rhea_reaction_side[-1]
            pattern1 = re.compile("R|L")
            if pattern1.match(rmc):
                gather_participant_data(description)

        elif "RHEA:" in node.text:
            first_entry_test = description.find("rh:substrates", ns) == None and description.find("rh:substratesOrProducts", ns) == None
            if first_entry_test:
                rhea_id = node.text
                rhea_entry = gather_rhea_data(rhea_id)
                rhea_lib.append(rhea_entry)

        else:
            continue
            
    for doc in rhea_lib:
        yield doc
