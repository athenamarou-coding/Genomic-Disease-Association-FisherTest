# mondo_utils.py
import json

MONDO_PREFIX = "http://purl.obolibrary.org/obo/MONDO_"


def part_1(path):
    with open(path, 'r') as f:
        data = json.load(f)
    graph = data.get("graphs", [{}])[0]
    nodes = graph.get("nodes", [])
    edges = graph.get("edges", [])

    if not nodes:
        print('No nodes found')
        return None

    disease_nodes = {}
    for n in nodes:
        node_id = n.get("id")
        if node_id and node_id.startswith(MONDO_PREFIX):
            short_id = node_id.split('/')[-1]
            disease_nodes[short_id] = {
                "id": short_id,
                "label": n.get("lbl"),
                "super_classes": set(),
                "sub_classes": set()
            }

    edge_count = 0
    for e in edges:
        pred = e.get("pred")

        if pred == "is_a":
            sub_url = e.get("sub")
            obj_url = e.get("obj")
            if (sub_url and sub_url.startswith(MONDO_PREFIX) and
                obj_url and obj_url.startswith(MONDO_PREFIX)):

                sub_id = sub_url.split('/')[-1]
                obj_id = obj_url.split('/')[-1]

                if sub_id in disease_nodes and obj_id in disease_nodes:
                    disease_nodes[sub_id]["super_classes"].add(obj_id)
                    disease_nodes[obj_id]["sub_classes"].add(sub_id)
                    edge_count += 1

    roots = []
    for node_id, node_data in disease_nodes.items():
        node_data["super_classes"] = list(node_data["super_classes"])
        node_data["sub_classes"] = list(node_data["sub_classes"])

        if not node_data["super_classes"]:
            roots.append(node_id)

    ontology_tree = {
        "nodes": disease_nodes,
        "roots": roots
    }

    return ontology_tree


def part_2(tree, start_mondo_id):
    nodes_dict = tree.get("nodes", {})

    if start_mondo_id not in nodes_dict:
        print(f"Σφάλμα: Το ID {start_mondo_id} δεν βρέθηκε στο δέντρο.")
        return []

    current_level = {start_mondo_id}
    visited = {start_mondo_id}
    depth = 0

    while depth < 3 and current_level:
        next_level = set()
        for nid in current_level:
            for parent_id in nodes_dict[nid].get("super_classes", []):
                if parent_id in nodes_dict and parent_id not in visited:
                    visited.add(parent_id)
                    next_level.add(parent_id)
        current_level = next_level
        depth += 1

    labels = []
    for nid in current_level:
        label = nodes_dict[nid].get("label")
        if label:
            labels.append(label)

    return sorted(labels)


def build_doid_to_mondo_mapping(path):
    """
    Γυρνάει το mondo.json και φτιάχνει:
      doid_to_mondo = { 'DOID:xxxx': ['MONDO_yyyy', ...], ... }
    """
    with open(path, "r", encoding="utf-8") as f:
        data = json.load(f)

    nodes = data["graphs"][0]["nodes"]

    doid_to_mondo = {}

    for n in nodes:
        node_id = n.get("id", "")
        if not node_id.startswith(MONDO_PREFIX):
            continue

        mondo_id = node_id.split("/")[-1]
        meta = n.get("meta", {})

        # 1. basicPropertyValues
        for bpv in meta.get("basicPropertyValues", []):
            val = bpv.get("val", "")
            if val.startswith("DOID:"):
                doid_to_mondo.setdefault(val, []).append(mondo_id)

        # 2. xrefs
        for xr in meta.get("xrefs", []):
            val = xr.get("val", "")
            if val.startswith("DOID:"):
                doid_to_mondo.setdefault(val, []).append(mondo_id)

    # Αφαίρεση διπλότυπων
    for doid, mondo_list in doid_to_mondo.items():
        doid_to_mondo[doid] = sorted(set(mondo_list))

    return doid_to_mondo
