'''PROJECT 2 - DISEASE ANALYSIS USING FISHER'S EXACT TEST AND MONDO ONTOLOGY'''

#Libraries
import csv
from collections import defaultdict
from scipy.stats import fisher_exact
#for bonus 2
from statsmodels.stats.multitest import multipletests 

#import from previous project
from mondo_utils import part_1, part_2, build_doid_to_mondo_mapping

#Paths
path_tsv = "/Users/athenamarounka/Documents/master/python/project_kantale/disease_model_annotations_fb_2025_02.tsv"
path_mondo = "/Users/athenamarounka/Documents/master/python/project_kantale/mondo.json"

#Dhmiourgia classes
class MondoKnowledge:
    '''Από το mondo ontology, φορτώνει το δέντρο και κάνει την αντιστοίχιση DOID σε Mondo Categories'''
    def __init__(self, json_path):
        self.json_path = json_path
        #Φορτώνω το δέντρο και το mapping από τα mondo_utils που έχω φτιάξει
        self.tree = part_1(self.json_path)
        self.doid_map = build_doid_to_mondo_mapping(self.json_path)
        
    def categories_for_doid(self, doid):
        '''Επιστρέφει ένα set με τις λατηγορίες Mondo για ένα DOID'''
        mondo_ids = self.doid_map.get(doid, [])
        categories = set()

        for mid in mondo_ids:
            #Με το part_2 από τα mondo_utils, βρίσκουμε κατηγορίες-προγόνους
            cats = part_2(self.tree, mid)
            categories.update(cats)
        return categories
    
    def get_mondo_ids(self, doid):
        return self.doid_map.get(doid, [])
    
class FlyBase:

    '''Ανάγνωση και φιλτράρισμα του FlyBase.tsv'''
    #Στήλες
    COL_FBGN_ID = 0
    COL_DO_QUALIFIER = 3
    COL_DO_ID = 4
    COL_DO_TERM = 5
    
    def __init__(self, tsv_path):
        self.tsv_path = tsv_path
        self.records = []
        
    def parse(self):
        '''Διαβάζει το αρχείο και βάζει αντικείμενα στην λίστα records'''
        results = []
        with open(self.tsv_path) as f:
            reader = csv.reader(f, delimiter= "\t")
            for row in reader:
                if not row:
                    continue
                if row[0].startswith('##'):
                    continue
                if len(row) <= self.COL_DO_TERM:
                    continue
                
                fbgn_id = row[self.COL_FBGN_ID].strip()
                do_qualifier = row[self.COL_DO_QUALIFIER].strip()
                do_id = row[self.COL_DO_ID].strip()
                do_term = row[self.COL_DO_TERM].strip()
                
                
                if not fbgn_id or not do_qualifier or not do_id:
                        continue

                results.append({
                    "fbgn_id": fbgn_id,
                    "do_qualifier": do_qualifier,
                    "do_id": do_id,
                    "do_term": do_term
                    })
        self.records = results
        print(f"{len(self.records)}")
        return self.records

class DiseaseAnalysis:
    """Στατιστική ανάλυση με δεδομένα από FlyBase and Mondo"""
    def __init__(self, parser: FlyBase, kb:MondoKnowledge):
        self.parser = parser
        self.kb = kb
        self.gene_info = {}
        self.all_qualifiers = set()
        self.all_categories = set()
        self.analysis_results = []


    def _prepare_data(self):
        '''Οργάνωση δεδομένων ανά γονίδιο'''
        for r in self.parser.records:
            fbgn = r["fbgn_id"]
            qualifier = r["do_qualifier"]
            doid = r["do_id"]
            
            #Ανάκτηση κατηγοριών από κλάση MondoKnowledge
            categories = self.kb.categories_for_doid(doid)
            if categories:
                if fbgn not in self.gene_info:
                    self.gene_info[fbgn] = {
                        "qualifiers": set(),
                        "mondo_categories": set()
                    }
                self.gene_info[fbgn]["qualifiers"].add(qualifier)
                self.gene_info[fbgn]["mondo_categories"].update(categories)
                self.all_qualifiers.add(qualifier)
                self.all_categories.update(categories)
                
        print(f"unique genes = {len(self.gene_info)}")
    
    def _compute_fisher(self, a, b, c, d):
        '''Υπολογισμός Fisher's Exact Test'''
        contingency_table = [[a, b], [c, d]]
        odds_ratio, p_value = fisher_exact(contingency_table, alternative='two-sided')
        return odds_ratio, p_value
    
    def run_analysis(self):
        '''Εκτέλεση ανάλυσης'''
        self._prepare_data()

        all_genes = list(self.gene_info.keys())
        sorted_categories = sorted(list(self.all_categories))
        sorted_qualifiers = sorted(list(self.all_qualifiers))
        total_genes = len(all_genes)
        
        for category in sorted_categories:
            for qualifier in sorted_qualifiers:
                
                #A: Has Category & Has Qualifier
                #B: Has Category & No Qualifier
                #C: No Category & Has Qualifier
                #D: No Category & No Qualifier

                count_A = 0
                count_B = 0
                count_C = 0
                count_D = 0
                
                for gene in all_genes:
                    g_data = self.gene_info[gene]
                    has_cat = category in g_data["mondo_categories"]
                    has_qual = qualifier in g_data["qualifiers"]
                    
                    if has_cat and has_qual:
                        count_A += 1
                    elif has_cat and not has_qual:
                        count_B += 1
                    elif not has_cat and has_qual:
                        count_C += 1
                    else:
                        count_D += 1

                if (count_A + count_B) > 0:
                    odds, p_val = self._compute_fisher(count_A, count_B, count_C, count_D)

                    #BONUS 1: EXPECTED VALUE AND FOLD CHANGE
                    expected_A = ((count_A + count_B) * (count_A + count_C)) / total_genes

                    #fold change = observed / expected
                    fold_change = (count_A / expected_A) if expected_A > 0 else float('inf')

                    self.analysis_results.append({
                        "p_value": p_val,
                        "odds_ratio": odds,
                        "category": category,
                        "qualifier": qualifier,
                        "table": [[count_A, count_B], [count_C, count_D]],

                        #metrics for bonus 1
                        "observed_A": count_A,
                        "expected_A": expected_A,
                        "fold_change": fold_change
                    })
    def print_top_results(self, top_n=10):
        '''Task 6: Εκτύπωση των top N αποτελεσμάτων με βάση το p-value'''
        print("\n" + "="*60)
        print(f"TASK 6 {top_n} Lowest P-values")
        print("="*60)
        
        sorted_results = sorted(self.analysis_results, key=lambda x: x["p_value"])
        
        for i, res in enumerate(sorted_results[:top_n]):
            print(f"{i+1}. Mondo Category: {res['category']}")
            print(f"   DO Qualifier:   {res['qualifier']}")
            print(f"   P-value:        {res['p_value']:.4e}")
            print(f"   Odds-ratio:     {res['odds_ratio']:.4f}")
            #BONUS 1 OUTPUT
            print(f"   Observed A:     {res['observed_A']}")
            print(f"   Expected A:     {res['expected_A']:.4f}")
            print(f"   Fold Change:    {res['fold_change']:.4f}")
            print("-" * 40)
            
    def run_multiple_testing_correction(self, alpha=0.05, method='fdr_bh'):
        '''Bonus 2: Διόρθωση πολλαπλών δοκιμών στα p-values'''
        
        p_values = [res['p_value'] for res in self.analysis_results]
        reject, pvals_corrected, _, _ = multipletests(p_values, alpha=alpha, method=method)
        
        for i, res in enumerate(self.analysis_results):
            res['p_value_corrected'] = pvals_corrected[i]
            res['reject_null'] = reject[i]
        
        print(f"\n --- BONUS 2 : Multiple Testing Correction ({method}) ---")
        
        rejected_count = 0
        
        for i, res in enumerate(self.analysis_results):
            if res['reject_null']:
                rejected_count += 1
                print(f"{i+1}. Mondo Category: {res['category']}")
                print(f"   DO Qualifier:   {res['qualifier']}")
                print(f"   Corrected P-value: {res['p_value_corrected']:.4e} (Reject Null Hypothesis)")
                print("-" * 40)
        print((f"Total Rejected Hypotheses: {rejected_count}"))

#EKTELESI
if __name__ == "__main__":
    mondo_kb = MondoKnowledge(path_mondo)
    flybase_parser = FlyBase(path_tsv)

    #Load FlyBase data
    flybase_parser.parse()

    #mondo_kb and flybase_parser inside analyzer
    analyzer = DiseaseAnalysis(flybase_parser, mondo_kb)
    analyzer.run_analysis()
    analyzer.print_top_results(10)
    analyzer.run_multiple_testing_correction()
