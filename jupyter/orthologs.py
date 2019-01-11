from SPARQLWrapper import SPARQLWrapper, JSON
from Bio.Graphics import GenomeDiagram
from Bio.Graphics.GenomeDiagram import CrossLink
from reportlab.lib import colors

from Bio.SeqFeature import SeqFeature, FeatureLocation
import pandas as pd
        
class SEARCH:

    def change_url(self, url):
       self.url = url['new']
       self.sparql = SPARQLWrapper(self.url)
       self.sparql.setReturnFormat(JSON)
        
    def __init__(self, url):
        #define url
        self.url = url
        #define sparql
        self.sparql = SPARQLWrapper(self.url)
        self.sparql.setReturnFormat(JSON)
        
    def url(self):
        print(self.url)
        
    def get_location(self, id):
        file = open("orthologs/gene_location.sparql", "r") 
        query = file.read()
        file.close()
        self.sparql.setQuery(query % id)
        # JSON example
        response = self.sparql.query().convert()
        result = []
        if response["results"]["bindings"]: 
            for item in response["results"]["bindings"]:
                result.append([
                item["gene_id"]["value"],
                item["location"]["value"],
                item["begin_ref"]["value"],
                item["begin_pos"]["value"],
                item["end_ref"]["value"],
                item["end_pos"]["value"]])
            df = pd.DataFrame(result)  
            df.columns = ["gene_id", "location", "begin_ref", "begin_pos", "end_ref", "end_pos" ]
            df = df.set_index("gene_id")
            df["begin_pos"] = pd.to_numeric(df["begin_pos"])
            df["end_pos"] = pd.to_numeric(df["end_pos"])
            return df 
        else:
            return pd.DataFrame()   
        
    def interval_genes(self, interval):
        file = open("orthologs/interval_genes.sparql", "r") 
        query = file.read()
        file.close()
        self.sparql.setQuery(query % {"beginRef" : interval.loc["begin"]["ref"], "beginPos" : interval.loc["begin"]["pos"], "endRef" : interval.loc["end"]["ref"], "endPos" : interval.loc["end"]["pos"]})
        # JSON example
        response = self.sparql.query().convert()
        result = []
        if response["results"]["bindings"]: 
            for item in response["results"]["bindings"]:
                row = []
                row.append(item["gene_id"]["value"])
                row.append(item["location"]["value"])
                row.append(item["begin_ref"]["value"])
                row.append(item["begin_pos"]["value"])
                row.append(item["end_ref"]["value"])
                row.append(item["end_pos"]["value"])
                if "ensembl_gene_id" in item.keys() :
                   row.append(item["ensembl_gene_id"]["value"])
                   row.append(item["ensembl_location"]["value"])
                   row.append(item["ensembl_begin_pos"]["value"])
                   row.append(item["ensembl_begin_ref"]["value"])
                   row.append(item["ensembl_end_pos"]["value"])
                   row.append(item["ensembl_end_ref"]["value"]) 
                   if "ortholog_gene_id" in item.keys() :
                     row.append(item["ortholog_gene_id"]["value"])
                     row.append(item["ortholog_location"]["value"])
                     row.append(item["ortholog_begin_pos"]["value"])
                     row.append(item["ortholog_begin_ref"]["value"])
                     row.append(item["ortholog_end_pos"]["value"])
                     row.append(item["ortholog_end_ref"]["value"])
                   else : 
                     row.append(None)
                     row.append(None)
                     row.append(None)
                     row.append(None)                
                     row.append(None)
                     row.append(None)
                else : 
                   row.append(None)
                   row.append(None)                
                   row.append(None)
                   row.append(None) 
                   row.append(None)
                   row.append(None)                
                   row.append(None)
                   row.append(None) 
                   row.append(None)
                   row.append(None)                
                   row.append(None)
                   row.append(None) 
                result.append(row)
            df = pd.DataFrame(result)  
            df.columns = ["gene_id", "location", "begin_ref", "begin_pos", "end_ref", "end_pos", "ensembl_gene_id", "ensembl_location", "ensembl_begin_pos", "ensembl_begin_ref", "ensembl_end_pos", "ensembl_end_ref", "ortholog_gene_id", "ortholog_location", "ortholog_begin_pos", "ortholog_begin_ref", "ortholog_end_pos", "ortholog_end_ref"]
            df["begin_pos"] = pd.to_numeric(df["begin_pos"])
            df["end_pos"] = pd.to_numeric(df["end_pos"])
            df["ensembl_begin_pos"] = pd.to_numeric(df["ensembl_begin_pos"])
            df["ensembl_end_pos"] = pd.to_numeric(df["ensembl_end_pos"])
            df["ortholog_begin_pos"] = pd.to_numeric(df["ortholog_begin_pos"])
            df["ortholog_end_pos"] = pd.to_numeric(df["ortholog_end_pos"])
            df = df.set_index("gene_id")
            return df 
        else:
            return pd.DataFrame()       
        
    def compute_interval(self, g1, g2):  
        locations = pd.concat([self.get_location(g1), self.get_location(g2)])
        display(locations[["location"]])
        if(len(locations.index)!=2) :
            print("unexpected number of rows in locations:",len(locations.index))
        elif(locations.iloc[0]['end_pos']>locations.iloc[1]['begin_pos']) :
            print("unexpected order",locations.index[0],"and",locations.index[1])
        else :
            result = []
            if locations.iloc[0]["end_pos"]>locations.iloc[0]["begin_pos"] :
              result.append(["begin", locations.iloc[0]["end_ref"], locations.iloc[0]["end_pos"]])
            else :
              result.append(["begin", locations.iloc[0]["begin_ref"], locations.iloc[0]["begin_pos"]])
            if locations.iloc[1]["begin_pos"]<locations.iloc[1]["end_pos"] :
              result.append(["end", locations.iloc[1]["begin_ref"], locations.iloc[1]["begin_pos"]])
            else :
              result.append(["end", locations.iloc[1]["end_ref"], locations.iloc[1]["end_pos"]])
            df = pd.DataFrame(result)
            df.columns = ["type", "ref", "pos" ]
            df = df.set_index("type")
            return df
    
    def diagram_genes(self, name, interval, genes):
        #initialise
        gd_diagram = GenomeDiagram.Diagram()        
        #add genes
        gd_track_for_features = gd_diagram.new_track(1,
                            name=name,
                            greytrack=True,
                            start=interval.loc["begin"]["pos"], end=interval.loc["end"]["pos"])
        gd_feature_set = gd_track_for_features.new_set()
        for gene in genes.groupby(["gene_id"]).first().itertuples():
            if gene.ensembl_gene_id is None:
              gene_seq_feature = SeqFeature(FeatureLocation(gene.begin_pos,gene.end_pos)) 
              gd_feature_set.add_feature(gene_seq_feature, color="lightblue", name=gene.Index, label=True)
            elif gene.ensembl_begin_pos > gene.ensembl_end_pos:
              gene_seq_feature = SeqFeature(FeatureLocation(int(gene.ensembl_end_pos),int(gene.ensembl_begin_pos)),strand=-1) 
              gd_feature_set.add_feature(gene_seq_feature, sigil="ARROW", color="blue", name=gene.Index, label=True)
            else:
              gene_seq_feature = SeqFeature(FeatureLocation(int(gene.ensembl_begin_pos),int(gene.ensembl_end_pos)),strand=+1) 
              gd_feature_set.add_feature(gene_seq_feature, sigil="ARROW", color="blue", name=gene.Index, label=True)
        #draw            
        gd_diagram.draw(format="linear", pagesize='A4', fragments=4,
                start=interval.loc["begin"]["pos"], end=interval.loc["end"]["pos"])        
        return gd_diagram
    
    def diagram_crosslinks(self, name, orthogonalName, interval, genes):
        #initialise
        gd_diagram = GenomeDiagram.Diagram()   
        minPos = interval.loc["begin"]["pos"]
        maxPos = interval.loc["end"]["pos"]  
        leftMargin = 10000
        rightMargin = 10000
        finalMaxPos = maxPos
        #add genes
        gd_track_for_genes = gd_diagram.new_track(1,
                            name=name,
                            greytrack=True,
                            start=minPos, end=maxPos)
        gd_feature_set = gd_track_for_genes.new_set()
        for gene in genes.groupby(["gene_id"]).first().itertuples():
            if gene.ensembl_gene_id is None:
              gene_seq_feature = SeqFeature(FeatureLocation(gene.begin_pos,gene.end_pos)) 
              gd_feature_set.add_feature(gene_seq_feature, color="lightblue", name=gene.Index, label=True)
            elif gene.ensembl_begin_pos > gene.ensembl_end_pos:
              gene_seq_feature = SeqFeature(FeatureLocation(int(gene.ensembl_end_pos),int(gene.ensembl_begin_pos)),strand=-1) 
              gd_feature_set.add_feature(gene_seq_feature, color="blue", name=gene.Index, label=True)
            else:
              gene_seq_feature = SeqFeature(FeatureLocation(int(gene.ensembl_begin_pos),int(gene.ensembl_end_pos)),strand=+1) 
              gd_feature_set.add_feature(gene_seq_feature, color="blue", name=gene.Index, label=True)
        #compute chromosomes for orthologs
        aggregations = {
            "ortholog_begin_pos" : ["min","max"],
            "ortholog_end_pos" : ["min","max"]
        }
        chromosomes = genes.groupby(["ortholog_begin_ref"]).agg(aggregations)
        for chromosome in chromosomes.itertuples():  
              minOrthologPos = min(chromosome._1,chromosome._2,chromosome._3,chromosome._4)
              maxOrthologPos = max(chromosome._1,chromosome._2,chromosome._3,chromosome._4)
              gd_track_for_orthologs = gd_diagram.new_track(1,
                            name=orthogonalName,
                            greytrack=True, scale_ticks=0, axis_labels=1,
                            start=minPos, end=min(maxPos, minPos+maxOrthologPos-minOrthologPos+leftMargin+rightMargin))                
              gd_feature_set = gd_track_for_orthologs.new_set() 
              positionShift = int(minPos - minOrthologPos + leftMargin)
              chromosomeGenes = genes.loc[genes["ortholog_begin_ref"] == chromosome.Index]
              for gene in chromosomeGenes.groupby(["ortholog_gene_id"]).first().itertuples():
                if gene.ortholog_begin_pos > gene.ortholog_end_pos:
                  gene_seq_feature = SeqFeature(FeatureLocation(int(gene.ortholog_end_pos)+positionShift,int(gene.ortholog_begin_pos)+positionShift),strand=-1) 
                  gd_feature_set.add_feature(gene_seq_feature, color="red", name=gene.Index, label=True)
                  link_xy = CrossLink((gd_track_for_genes, int(gene.begin_pos), int(gene.end_pos)),
                            (gd_track_for_orthologs, int(gene.ortholog_end_pos)+positionShift,int(gene.ortholog_begin_pos)+positionShift),
                            colors.lightgrey, colors.grey)
                  gd_diagram.cross_track_links.append(link_xy)  
                else:
                  gene_seq_feature = SeqFeature(FeatureLocation(int(gene.ortholog_begin_pos)+positionShift,int(gene.ortholog_end_pos)+positionShift),strand=+1) 
                  gd_feature_set.add_feature(gene_seq_feature, color="red", name=gene.Index, label=True)
                  link_xy = CrossLink((gd_track_for_genes, int(gene.begin_pos), int(gene.end_pos)),
                            (gd_track_for_orthologs, int(gene.ortholog_begin_pos)+positionShift,int(gene.ortholog_end_pos)+positionShift),
                            colors.lightgrey, colors.grey)
                  gd_diagram.cross_track_links.append(link_xy)  
              for gene in chromosomeGenes.groupby(["gene_id"]).first().itertuples():  
                if not (gene.ortholog_gene_id is None):                  
                    if gene.ortholog_begin_pos > gene.ortholog_end_pos:
                      link_xy = CrossLink((gd_track_for_genes, int(gene.begin_pos), int(gene.end_pos)),
                                (gd_track_for_orthologs, int(gene.ortholog_end_pos)+positionShift,int(gene.ortholog_begin_pos)+positionShift),
                                colors.lightgrey, colors.grey)
                      gd_diagram.cross_track_links.append(link_xy)  
                    else:
                      link_xy = CrossLink((gd_track_for_genes, int(gene.begin_pos), int(gene.end_pos)),
                                (gd_track_for_orthologs, int(gene.ortholog_begin_pos)+positionShift,int(gene.ortholog_end_pos)+positionShift),
                                colors.lightgrey, colors.grey)
                      gd_diagram.cross_track_links.append(link_xy) 
        
        #draw            
        gd_diagram.draw(format="linear", pagesize='A4', fragments=1,
                start=minPos, end=finalMaxPos)        
        return gd_diagram
     
        

    
    

    