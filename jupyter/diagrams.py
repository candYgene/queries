from Bio.Graphics import GenomeDiagram
from Bio.Graphics.GenomeDiagram import CrossLink
from reportlab.lib import colors

from Bio.SeqFeature import SeqFeature, FeatureLocation

class DIAGRAMS:
    
    def crosslinks(name, orthogonalName, interval, genes):
        #initialise
        gd_diagram = GenomeDiagram.Diagram()   
        minPos = interval["begin_pos"]
        maxPos = interval["end_pos"]  
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
            if gene.begin_pos > gene.end_pos:
                gene_seq_feature = SeqFeature(FeatureLocation(int(gene.end_pos),int(gene.begin_pos)),strand=-1) 
                gd_feature_set.add_feature(gene_seq_feature, color="blue", name=gene.Index, label=True)
            else:
                gene_seq_feature = SeqFeature(FeatureLocation(int(gene.begin_pos),int(gene.end_pos)),strand=+1) 
                gd_feature_set.add_feature(gene_seq_feature, color="blue", name=gene.Index, label=True)
        #compute chromosomes for orthologs
        aggregations = {
            "ortholog_begin_pos" : ["min","max"],
            "ortholog_end_pos" : ["min","max"]
        }
        chromosomes = genes.groupby(["ortholog_taxon_id", "ortholog_chrom"]).agg(aggregations)
        for chromosome in chromosomes.itertuples():  
              minOrthologPos = min(chromosome._1,chromosome._2,chromosome._3,chromosome._4)
              maxOrthologPos = max(chromosome._1,chromosome._2,chromosome._3,chromosome._4)
              gd_track_for_orthologs = gd_diagram.new_track(1,
                            name=orthogonalName,
                            greytrack=True, scale_ticks=0, axis_labels=1,
                            start=minPos, end=min(maxPos, minPos+maxOrthologPos-minOrthologPos+leftMargin+rightMargin))                
              gd_feature_set = gd_track_for_orthologs.new_set() 
              positionShift = int(minPos - minOrthologPos + leftMargin)
              chromosomeGenes = genes.loc[(genes["ortholog_taxon_id"] == chromosome.Index[0]) & 
                                          (genes["ortholog_chrom"] == chromosome.Index[1])]
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
         