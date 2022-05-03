import re
import xml.etree.ElementTree as ET
from xml.etree.ElementTree import Element
from typing import List


# dictionary with namespaces
NAMESPACES = {
    "rh": "http://rdf.rhea-db.org/",
    "rdfs": "http://www.w3.org/2000/01/rdf-schema#",
    "rdf": "http://www.w3.org/1999/02/22-rdf-syntax-ns#"
}

# Collection of commonly used attribute keys to `xml.etree.ElementTree.Element.attrib`
ATTRIB_KEYS = {
    # expands to "{http://www.w3.org/1999/02/22-rdf-syntax-ns#}resource"
    "rdf:resource": f"{{{NAMESPACES['rdf']}}}resource",
    # expands to "{http://www.w3.org/1999/02/22-rdf-syntax-ns#}about"
    "rdf:about": f"{{{NAMESPACES['rdf']}}}about"
}

# Collection of commonly used base IRIs
BASE_IRIS = {
    "rh": NAMESPACES["rh"],
    "pubmed": "http://rdf.ncbi.nlm.nih.gov/pubmed/",
    "obo": "http://purl.obolibrary.org/obo/",
    "ec": "http://purl.uniprot.org/enzyme/"
}


class ReactivePartDataFactory:
    """
    A macromolecule (of type "rh:GenericCompound") has reactive parts, specified by description like:

        <rdf:Description rdf:about="http://rdf.rhea-db.org/Compound_9846_rp1">
            <rdfs:subClassOf rdf:resource="http://rdf.rhea-db.org/ReactivePart"/>
            <rh:name>N(6)-methyl-L-lysine residue</rh:name>
            <rh:htmlName>...</rh:htmlName>
            <rh:formula>C7H15N2O</rh:formula>
            <rh:charge>1</rh:charge>
            <rh:chebi rdf:resource="http://purl.obolibrary.org/obo/CHEBI_61929"/>
            <rdfs:subClassOf rdf:resource="http://purl.obolibrary.org/obo/CHEBI_61929"/>
        </rdf:Description>

    The relative IRI "Compound_9846_rp1" indicates this is the first reactive part of compound "9846".

    This data factory produce all such associations in tuples (comp_num, rp_entry).
    """
    relative_iri_pattern = re.compile(r"Compound_\d+_rp\d")

    @classmethod
    def is_valid_relative_iri(cls, relative_iri: str):
        match = cls.relative_iri_pattern.fullmatch(relative_iri)
        return match is not None

    # Adds name key and associated value to compounds with reactive parts
    @classmethod
    def _add_rp_name(cls, rp_entry: dict, description: Element):
        node = description.find("rh:name", NAMESPACES)
        if node is not None:
            rp_name = node.text
            rp_entry["name"] = rp_name

    # Adds formula key and associated value to compounds with reactive parts
    @classmethod
    def _add_rp_formula(cls, rp_entry: dict, description: Element):
        node = description.find("rh:formula", NAMESPACES)
        if node is not None:
            rp_formula = node.text
            rp_entry["formula"] = rp_formula

    # Adds charge key and associated value to compounds with reactive parts
    @classmethod
    def _add_rp_charge(cls, rp_entry: dict, description: Element):
        node = description.find("rh:charge", NAMESPACES)
        if node is not None:
            rp_charge = int(node.text)
            rp_entry["charge"] = rp_charge

    # Adds chebi key and associated id to compounds with reactive parts
    @classmethod
    def _add_rp_chebi_id(cls, rp_entry: dict, description: Element):
        node = description.find("rh:chebi", NAMESPACES)
        if node is not None:
            rp_chebi = node.attrib[ATTRIB_KEYS["rdf:resource"]].lstrip(BASE_IRIS["obo"]).replace("_", ":")
            rp_entry["chebi_id"] = rp_chebi

    # Adds reactive part information to the list associated with the key "reactive_parts" to compounds with this annoted
    @classmethod
    def produce(cls, relative_iri: str, description: Element):
        # we can assume here that relative IRI has a pattern of "Compound_\d+_rp\d", e.g. "Compound_10594_rp2"
        # comp_num = re.sub(pattern=r"_rp\d", repl="", string=relative_iri).lstrip("Compound_")
        comp_num = relative_iri.split("_")[1]

        rp_entry = {}
        cls._add_rp_name(rp_entry, description)
        cls._add_rp_formula(rp_entry, description)
        cls._add_rp_charge(rp_entry, description)
        cls._add_rp_chebi_id(rp_entry, description)

        yield comp_num, rp_entry


class SideDataFactory:
    """
    Each reaction specifies its two side ("L" for left, "R" for right). E.g.

        <rdf:Description rdf:about="http://rdf.rhea-db.org/35975">
            # some other tags ignored
            <rh:side rdf:resource="http://rdf.rhea-db.org/35975_L"/>
        </rdf:Description>

        <rdf:Description rdf:about="http://rdf.rhea-db.org/35975">
            <rh:side rdf:resource="http://rdf.rhea-db.org/35975_R"/>
        </rdf:Description>

    Each side description specifies a "contains" relationship for each participating compound. E.g.

        <rdf:Description rdf:about="http://rdf.rhea-db.org/35975_R">
            <rh:contains rdf:resource="http://rdf.rhea-db.org/Participant_35975_compound_3512"/>
            <rh:contains1 rdf:resource="http://rdf.rhea-db.org/Participant_35975_compound_3512"/>
        </rdf:Description>

    Note that "rh:contains" only specifies containing without stoichiometry. "rh:contains1" indicates a containing
    relationship along with stoichiometry 1.

    Therefore from this side description we derive the following associations:

    - The right side of reaction "35975" contains compound "3512", and
    - the compound's stoichiometry is 1 in the reaction.

    We can also get the location ("in" or "out") of a compound, if its reaction is a transport reaction and the
    compound is found in both sides of the reaction. E.g.

        <rdf:Description rdf:about="http://rdf.rhea-db.org/69560_L">
            <rh:contains rdf:resource="http://rdf.rhea-db.org/Participant_69560_compound_3249_out"/>
            <rh:contains1 rdf:resource="http://rdf.rhea-db.org/Participant_69560_compound_3249_out"/>
        </rdf:Description>

    The location of compound "3249" is "out".

    This data factory produces association quadruples of (rhea_id, side_key, comp_num, participation_entry), which
    indicates:

    - The reaction with `rhea_id`, on its left/right side specified by `side_key`, contains a compound indicated by
    `comp_num`, and
    - the compound's stoichiometry and location (if any) are wrapped in the `participation_entry`.
    """

    # expands to "{http://rdf.rhea-db.org/}contains"
    contains_prefix = f"{{{NAMESPACES['rh']}}}contains"

    side_char_to_key = {
        "L": "side_l",
        "R": "side_r"
    }

    @classmethod
    def is_valid_relative_iri(cls, relative_iri: str):
        rmc = relative_iri[-1]  # right-most character
        return rmc in cls.side_char_to_key

    @classmethod
    def is_stoichiometric_tag(cls, tag: str):
        """
        Check if a tag has the form of "rh:contains[x]" from which the stoichiometry of [x] can be inferred
        """
        return (tag != cls.contains_prefix) and tag.startswith(cls.contains_prefix)

    @classmethod
    def get_stoichiometry(cls, tag: str):
        """
        Trim the prefix "rh:contains" from a stoichimetric tag to get the stoichiometry.

        Note that there exist special stoichimetric tags like

        - "containsN",
        - "contains2n" (no idea why it's not "2N"),
        - "containsNplus1", and
        - "containsNminus1"

        Therefore stoichiometry's datatype must be string.
        """
        return tag.lstrip(cls.contains_prefix)

    @classmethod
    def produce(cls, relative_iri: str, description: Element):
        # we assume relative_iri is valid
        side_char = relative_iri[-1]
        side_key = cls.side_char_to_key[side_char]

        child_tags = [child.tag for child in description if cls.is_stoichiometric_tag(child.tag)]
        for tag in child_tags:
            # namespaces not needed here for description.find() because `tag` has the expanded namespace already
            contained_absolute_iri = description.find(tag).attrib[ATTRIB_KEYS["rdf:resource"]]
            contained_relative_iri = contained_absolute_iri.lstrip(BASE_IRIS['rh'])

            # E.g. contained_relative_iri = "Participant_69560_compound_3249_out"
            contained_relative_iri_parts = contained_relative_iri.split("_")
            rhea_num = contained_relative_iri_parts[1]
            comp_num = contained_relative_iri_parts[3]
            location = contained_relative_iri_parts[4] if len(contained_relative_iri_parts) == 5 else None

            """
            Here we assume that rhea id can be inferred from the contained IRI (because accession id is not
                available).
            E.g. an IRI like "Participant_35975_compound_3512" indicates that compound_lib["3512"] is associated
                with reaction "RHEA:35975"
            It's also possible to infer from the upper-level IRI (the `relative_iri` argument). E.g. "35975_R" in
                the above example
            """
            rhea_id = "RHEA:" + rhea_num
            stoich = cls.get_stoichiometry(tag)

            participation_entry = {"stoich": stoich}

            # Adds positional information for compounds to rhea entries that specify a transport reaction
            if location:
                participation_entry["location"] = location

            yield rhea_id, side_key, comp_num, participation_entry


class CompoundDataFactory:
    """
    A compound's description is like

        <rdf:Description rdf:about="http://rdf.rhea-db.org/Compound_1454">
            <rh:id rdf:datatype="http://www.w3.org/2001/XMLSchema#long">1454</rh:id>
            <rh:accession>CHEBI:58413</rh:accession>
            <rh:name>(R)-6-hydroxynicotine</rh:name>
            <rh:htmlName>(&lt;i&gt;R&lt;/i&gt;)-6-hydroxynicotine</rh:htmlName>
            <rh:formula>C10H15N2O</rh:formula>
            <rh:charge rdf:datatype="http://www.w3.org/2001/XMLSchema#int">1</rh:charge>
            <rdfs:subClassOf rdf:resource="http://rdf.rhea-db.org/SmallMolecule"/>
            <rh:chebi rdf:resource="http://purl.obolibrary.org/obo/CHEBI_58413"/>
            <rdfs:subClassOf rdf:resource="http://purl.obolibrary.org/obo/CHEBI_58413"/>
        </rdf:Description>

    Note that only such descriptions with a valid accession ID will be parsed by this data factory. Currently there
    are 3 types of valid accession IDs for compounds, i.e. "CHEBI", "GENERIC", and "POLYMER".

    This data factory produces basic information (without reactive parts, stoichiometry, nor location) for each
    compound.
    """
    compound_prefixes = set(["CHEBI:", "GENERIC:", "POLYMER:"])

    @classmethod
    def is_valid_accession_id(cls, accession_id: str):
        for prefix in cls.compound_prefixes:
            if accession_id.startswith(prefix):
                return True
        return False

    # Adds ID key and value to compound entries
    @classmethod
    def _add_comp_id(cls, comp_entry: dict, accession_id: str):
        if "CHEBI:" in accession_id:
            comp_entry["chebi_id"] = accession_id
        elif "GENERIC:" in accession_id:
            comp_entry["generic_id"] = accession_id.lstrip("GENERIC:")
        elif "POLYMER:" in accession_id:
            comp_entry["poly_id"] = accession_id.lstrip("POLYMER:")
        else:
            raise ValueError(f"Cannot recognize accession type. Got accession id {accession_id}")

    # Adds name key and id to compound entries
    @classmethod
    def _add_comp_name(cls, comp_entry: dict, description: Element):
        node = description.find("rh:name", NAMESPACES)
        if node is not None:
            comp_name = node.text
            comp_entry["name"] = comp_name

    # Adds formula key and value to compound entries
    @classmethod
    def _add_comp_formula(cls, comp_entry: dict, description: Element):
        node = description.find("rh:formula", NAMESPACES)
        if node is not None:
            formula = node.text
            if formula is not None:
                formula = formula.rstrip("<i><sub>n</sub></i>")
                comp_entry["formula"] = formula

    # Adds charge key and value to compound entries
    @classmethod
    def _add_comp_charge(cls, comp_entry: dict, description: Element):
        node = description.find("rh:charge", NAMESPACES)
        if node is not None:
            comp_charge = node.text.rstrip("<i><sub>n</sub></i>")
            # comp_charge can be a string like '(-4)(-1)' so its datatype cannot be integer
            comp_entry["charge"] = comp_charge

    @classmethod
    def produce(cls, relative_iri: str, accession_id: str, description: Element):
        comp_entry = {}

        # we can assume here that relative IRI has a pattern of "Compound_\d", e.g. "Compound_10594"
        comp_num = relative_iri.split("_")[1]
        comp_entry["comp_num"] = comp_num

        cls._add_comp_id(comp_entry, accession_id)
        cls._add_comp_name(comp_entry, description)
        cls._add_comp_formula(comp_entry, description)
        cls._add_comp_charge(comp_entry, description)

        yield comp_entry

    @classmethod
    def pack(cls, comp_entries: List[dict]):
        """
        Pack a list of compound entries into a dictionary of {comp_num : comp_entry}
        """
        return dict((comp_entry["comp_num"], comp_entry) for comp_entry in comp_entries)


class ReactionDataFactory:
    """
    A reaction's description is like

        <rdf:Description rdf:about="http://rdf.rhea-db.org/10000">
            <rdfs:subClassOf rdf:resource="http://rdf.rhea-db.org/Reaction"/>
            <rh:id rdf:datatype="http://www.w3.org/2001/XMLSchema#long">10000</rh:id>
            <rh:accession>RHEA:10000</rh:accession>
            <rdfs:label>H2O + pentanamide = NH4(+) + pentanoate</rdfs:label>
            <rh:equation>H2O + pentanamide = NH4(+) + pentanoate</rh:equation>
            <rh:htmlEquation>...</rh:htmlEquation>
            <rh:directionalReaction rdf:resource="http://rdf.rhea-db.org/10001"/>
            <rh:directionalReaction rdf:resource="http://rdf.rhea-db.org/10002"/>
            <rh:bidirectionalReaction rdf:resource="http://rdf.rhea-db.org/10003"/>
            <rh:status rdf:resource="http://rdf.rhea-db.org/Approved"/>
            <rh:isChemicallyBalanced rdf:datatype="http://www.w3.org/2001/XMLSchema#boolean">
                true
            </rh:isChemicallyBalanced>
            <rh:isTransport rdf:datatype="http://www.w3.org/2001/XMLSchema#boolean">false</rh:isTransport>
            <rdfs:comment>...</rdfs:comment>
            <rh:ec rdf:resource="http://purl.uniprot.org/enzyme/3.5.1.50"/>
            <rdfs:seeAlso rdf:resource="http://identifiers.org/biocyc/METACYC:PENTANAMIDASE-RXN"/>
            <rdfs:seeAlso rdf:resource="http://purl.obolibrary.org/obo/GO_0050168"/>
            <rh:side rdf:resource="http://rdf.rhea-db.org/10000_L"/>
        </rdf:Description>

    Note that only such descriptions with a valid accession ID will be selected to this data factory. Currently there
    are 1 type of valid accession IDs for reactions, i.e. those starting with "RHEA".

    Also note that for each reaction, there will be 4 variants, i.e.

    - the master reaction (direction undefined, e.g. "RHEA:10000")
    - 2 directional reactions (left-to-right, e.g. "RHEA:10001", and right-to-left, e.g. "RHEA:10002")
    - the bidirectional reaction (e.g. RHEA:10003)

    In this data factory, only the master reactions will be parsed to individual entries, the directional and 
    bidirectional reactions will be attached to their master reactions as "children_rheas".
    """
    @classmethod
    def is_valid_accession_id(cls, accession_id):
        return accession_id.startswith("RHEA:")

    @classmethod
    def is_master_reaction(cls, description):
        """
        There are multiple ways to tell if a reaction is a master reaction.

        Method 1: tell by "rdfs:subClassOf", whose values are "Reaction", "DirectionalReaction",
            "BidirectionalReaction".
        Method 2: tell by the existence of "rh:substrates" and/or "rh:products" (only in "DirectionalReaction"),
            plus "rh:substratesOrProducts" (only in "BidirectionalReaction")

        Here we use method 2.

        Note that a master reaction's RHEA ID is not necessarily a multiple of 4. E.g. "RHEA:26018" for some reason is
        not used, and the next master reaction is "RHEA:26019". Therefore there is no modulo relationship between a
        reaction's type and its RHEA ID.
        """
        return (description.find("rh:substrates", NAMESPACES) is None) and \
            (description.find("rh:substratesOrProducts", NAMESPACES) is None)

    # Adds equation key and associated value to reaction entry
    @classmethod
    def _add_rhea_equation(cls, reaction_entry: dict, description: Element):
        node = description.find("rh:equation", NAMESPACES)
        if node is not None:
            reaction_entry["equation"] = node.text

    # Adds is_transport key and associated boolean to reaction entry
    @classmethod
    def _add_rhea_transport(cls, reaction_entry: dict, description: Element):
        node = description.find("rh:isTransport", NAMESPACES)
        if node is not None:
            """
            "rh:isTransport" has specified rdf:datatype="http://www.w3.org/2001/XMLSchema#boolean",
            therefore only 2 unique values are possible, "true" and "false"
            """
            is_transport = node.text  # string type
            is_transport = (is_transport == "true")  # boolean type
            reaction_entry["is_transport"] = is_transport

    # Adds ec_link and ec_id keys and associated values to reaction entry
    # ENZYME is an enzyme nomenclature database, which assigns an EC (Enzyme Commission) number for each enzyme
    @classmethod
    def _add_rhea_ec(cls, reaction_entry: dict, description: Element):
        node = description.find("rh:ec", NAMESPACES)
        if node is not None:
            ec_link = node.attrib[ATTRIB_KEYS["rdf:resource"]]
            ec_id = ec_link.lstrip(BASE_IRIS["ec"])

            reaction_entry["ec_link"] = ec_link
            reaction_entry["ec_id"] = ec_id

    # Adds status key and associated value to reaction entry. 3 possible values: Approved, Preliminary, Obsolete
    @classmethod
    def _add_rhea_status(cls, reaction_entry: dict, description: Element):
        node = description.find("rh:status", NAMESPACES)
        if node is not None:
            status = node.attrib[ATTRIB_KEYS["rdf:resource"]].lstrip(BASE_IRIS["rh"])
            reaction_entry["status"] = status

    # Adds citations key and associated values in a list to reaction entry.
    #  Some entries will have no citations and thus no citations key
    @classmethod
    def _add_rhea_citations(cls, reaction_entry: dict, description: Element):
        nodes = description.findall("rh:citation", NAMESPACES)
        if nodes:
            for node in nodes:
                citation = node.attrib[ATTRIB_KEYS["rdf:resource"]].lstrip(BASE_IRIS["pubmed"])
                citation = "PMID:" + citation
                reaction_entry.setdefault("citations", []).append(citation)

    # Adds the children_rheas key and list of associated rhea ids (should be 3) to reaction entry
    @classmethod
    def _add_rhea_children(cls, reaction_entry: dict, description: Element):
        directional_reactions = description.findall("rh:directionalReaction", NAMESPACES)
        if directional_reactions:
            for reaction in directional_reactions:
                child_absoulte_iri = reaction.attrib[ATTRIB_KEYS["rdf:resource"]]
                child_relative_iri = child_absoulte_iri.lstrip(BASE_IRIS["rh"])
                child_rhea_id = "RHEA:" + child_relative_iri

                reaction_entry.setdefault("children_rheas", []).append(child_rhea_id)

        bidirectional_reaction = description.find("rh:bidirectionalReaction", NAMESPACES)
        if bidirectional_reaction is not None:
            child_absoulte_iri = bidirectional_reaction.attrib[ATTRIB_KEYS["rdf:resource"]]
            child_relative_iri = child_absoulte_iri.lstrip(BASE_IRIS["rh"])
            child_rhea_id = "RHEA:" + child_relative_iri

            reaction_entry.setdefault("children_rheas", []).append(child_rhea_id)

    # Fills rhea entries with associated information
    @classmethod
    def produce(cls, accession_id: str, description: Element):
        reaction_entry = {}

        reaction_entry["rhea_id"] = accession_id

        # reaction_entry["side_l"] = []
        # reaction_entry["side_r"] = []

        cls._add_rhea_equation(reaction_entry, description)
        cls._add_rhea_transport(reaction_entry, description)
        cls._add_rhea_ec(reaction_entry, description)
        cls._add_rhea_status(reaction_entry, description)
        cls._add_rhea_citations(reaction_entry, description)
        cls._add_rhea_children(reaction_entry, description)

        yield reaction_entry

    @classmethod
    def pack(cls, reaction_entries: List[dict]):
        return dict((reaction_entry["rhea_id"], reaction_entry) for reaction_entry in reaction_entries)


def load_annotations(data_folder):
    """
    Using ElementTree, rhea.rdf is collapsed into a hierarchy of tags with associated information accessed with .find()
    and .findall() functions.

    The main for-loop catches reactive part associations, reaction side associations, compound entries, and reaction
    entries into 4 lists.

    Reactive parts are augmented to their associated compounds in the 2nd for-loop.

    The 3rd for-loop creates side components from compounds and participation fields, and then attaches side components
    to the associated reactions.

    The 4th for-loop yields all reaction docments.
    """
    rhea_rdf = open(data_folder + "/rhea.rdf", "r")
    tree = ET.parse(rhea_rdf)
    root = tree.getroot()

    reactive_part_associations = []
    side_associations = []
    compound_entries = []
    reaction_entries = []

    for description in root.findall("rdf:Description", NAMESPACES):
        absolute_iri = description.attrib[ATTRIB_KEYS["rdf:about"]]
        relative_iri = absolute_iri.lstrip(BASE_IRIS["rh"])

        accession = description.find("rh:accession", NAMESPACES)

        if accession is None:
            if ReactivePartDataFactory.is_valid_relative_iri(relative_iri=relative_iri):
                for rp_assoc in ReactivePartDataFactory.produce(relative_iri=relative_iri, description=description):
                    reactive_part_associations.append(rp_assoc)
            elif SideDataFactory.is_valid_relative_iri(relative_iri=relative_iri):
                for side_assoc in SideDataFactory.produce(relative_iri=relative_iri, description=description):
                    side_associations.append(side_assoc)
        else:
            accession_id = accession.text
            if CompoundDataFactory.is_valid_accession_id(accession_id=accession_id):
                for comp_entry in CompoundDataFactory.produce(relative_iri=relative_iri,
                                                              accession_id=accession_id,
                                                              description=description):
                    compound_entries.append(comp_entry)
            elif ReactionDataFactory.is_valid_accession_id(accession_id=accession_id) and \
                    ReactionDataFactory.is_master_reaction(description=description):
                for rhea_entry in ReactionDataFactory.produce(accession_id=accession_id, description=description):
                    reaction_entries.append(rhea_entry)

    compound_lib = CompoundDataFactory.pack(comp_entries=compound_entries)
    reaction_lib = ReactionDataFactory.pack(reaction_entries=reaction_entries)

    for comp_num, rp_entry in reactive_part_associations:
        compound_lib[comp_num].setdefault("reactive_parts", []).append(rp_entry)

    for rhea_id, side_key, comp_num, participation_entry in side_associations:
        comp_entry = compound_lib[comp_num]
        if "comp_num" in comp_entry:
            del comp_entry["comp_num"]

        side_component = {
            **comp_entry,  # the "participant"
            **participation_entry  # describes how this compound partipates in the reaction
        }
        reaction_lib[rhea_id].setdefault(side_key, []).append(side_component)

    for reaction_entry in reaction_lib.values():
        reaction_entry["_id"] = reaction_entry["rhea_id"]
        del reaction_entry["rhea_id"]
        yield reaction_entry
