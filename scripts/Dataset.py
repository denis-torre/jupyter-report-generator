#######################################################
########## 1. Libraries ###############################
#######################################################
import gzip, urllib2, os, requests, json, qgrid
import pandas as pd
import numpy as np
from StringIO import StringIO
from clustergrammer import Network
from sklearn.decomposition import PCA
import plotly.graph_objs as go
from plotly.offline import iplot
from clustergrammer_widget import *
from time import sleep
from IPython.display import HTML

from rpy2.robjects import r, pandas2ri
r.source('../scripts/Dataset.R')
pandas2ri.activate()
pd.set_option('max.colwidth', -1)

#################################################################
#################################################################
########## 1. Dataset Class #####################################
#################################################################
#################################################################

#######################################################
########## 1. Initialize ##############################
#######################################################

class Dataset:

    def __init__(self, matrix_url):

        # Print
        print 'Fetching expression data from ARCHS4...'

        # Open URL
        request = urllib2.Request(matrix_url)

        # Add header
        request.add_header('Accept-encoding', 'gzip')

        # Read response
        response = urllib2.urlopen(request)

        # Convert response
        buf = StringIO(response.read())

        # Open gzip file
        f = gzip.GzipFile(fileobj=buf)

        # Read data
        data = f.read()

        # Get platform
        self.platform_accession = os.path.basename(matrix_url).split('_')[1]

        # Check platform
        if self.platform_accession in ['GPL11154', 'GPL13112', 'GPL16791', 'GPL17021']:

            # Get expression dataframe
            rawcount_dataframe = pd.DataFrame([x.split('\t') for x in data.split('\n')[1:] if '!' not in x])

            # Fix axis names
            self.rawcount_dataframe = rawcount_dataframe.rename(columns=rawcount_dataframe.iloc[0]).drop(0).set_index('ID_REF').fillna(0).astype('int')

            # Get metadata dataframe
            metadata_dataframe = pd.DataFrame([x.split('\t') for x in data.split('\n')[1:] if any(y in x for y in ['!Sample', '!^SAMPLE'])]).set_index(0)

            # Get title conversion
            self.sample_title_dict = {sample_accession: '{sample_title} ({sample_accession})'.format(**locals()) for sample_accession, sample_title in metadata_dataframe.loc[['!^SAMPLE', '!Sample_title']].T.as_matrix() if sample_accession}

            # Get metadata dict
            metadata_dict = [{term_string.split(': ')[0]: term_string.split(': ')[1] for term_string in term_list} for term_list in np.array_split(metadata_dataframe.loc['!Sample_characteristics_ch1'].dropna().tolist(), len(self.sample_title_dict.keys()))]
            
            # Create dict
            self.sample_metadata_dataframe = pd.DataFrame({sample_accession: metadata_dict for sample_accession, metadata_dict in zip(metadata_dataframe.loc['!^SAMPLE'], metadata_dict)}).T

        # Print
        print 'Done!'

#######################################################
########## 2. Plot PCA ################################
#######################################################

    def plot_pca(self):

        # Perform PCA
        pca = PCA(n_components=3)
        pca.fit(np.log10(self.rawcount_dataframe+1))

        # Plot it
        trace = go.Scatter3d(
            x=pca.components_[0],
            y=pca.components_[1],
            z=pca.components_[2],
            mode='markers',
            hoverinfo='text',
            text=[self.sample_title_dict[x] for x in self.rawcount_dataframe.columns],
            marker=dict(
                size=12,
                colorscale='Viridis',   # choose a colorscale
                opacity=0.8
            )
        )

        data = [trace]
        go.layout = go.Layout(
            title='PCA',
            margin=dict(
                l=0,
                r=0,
                b=0,
                t=0
            ),
            scene = dict(
                xaxis = dict(
                    title = 'PC1 ({:0.2f}% var.)'.format(pca.explained_variance_ratio_[0]*100)
                ),
                yaxis = dict(
                    title = 'PC2 ({:0.2f}% var.)'.format(pca.explained_variance_ratio_[1]*100)
                ),
                zaxis = dict(
                    title = 'PC3 ({:0.2f}% var.)'.format(pca.explained_variance_ratio_[2]*100)
                )
            )
        )
        fig = go.Figure(data=data, layout=go.layout)
        return iplot(fig)

#######################################################
########## 3. Plot Clustergram ########################
#######################################################

    def plot_clustergram(self, n_genes=500):

        # Get top variable genes
        top_variance_genes = self.rawcount_dataframe.apply(np.var, 1).nlargest(n_genes).index.tolist()

        # Rename
        clustergrammer_dataframe = self.rawcount_dataframe.rename(columns=self.sample_title_dict).loc[top_variance_genes]

        # Z-score
        clustergrammer_dataframe = ((clustergrammer_dataframe.T - clustergrammer_dataframe.T.mean())/clustergrammer_dataframe.T.std()).T

        # Add cats
        sample_cats = [{'title': index, 'cats':{value:[self.sample_title_dict[x] for x in rowData[rowData==value].index.tolist()] for value in set(rowData.dropna().values)}} for index, rowData in self.sample_metadata_dataframe.T.iterrows()]
        
        # Create network
        net = Network(clustergrammer_widget)

        # Load file
        net.load_df(clustergrammer_dataframe)

        # Add categories
        try:
            net.add_cats('col', sample_cats)
        except:
            pass

        # Cluster
        net.cluster()

        # Return widget
        return net.widget()

#######################################################
########## 4. Set Experimental Design #################
#######################################################

    def set_experimental_design(self, experimental_samples, control_samples):

        # Add samples
        sample_dict = {'experimental': experimental_samples, 'control': control_samples}

        # Create dataframe
        self.design_dataframe = pd.DataFrame({group_label: {sample:int(sample in group_samples) for sample in self.rawcount_dataframe.columns} for group_label, group_samples in sample_dict.iteritems()})

#######################################################
########## 5. Run limma ###############################
#######################################################

    def run_limma(self):

        # Print
        print 'Running differential expression analysis...'

        # Run limma
        self.limma_dataframe = pandas2ri.ri2py(r.run_limma(pandas2ri.py2ri(self.rawcount_dataframe), pandas2ri.py2ri(self.design_dataframe)))

        # Filter
        self.limma_dataframe = self.limma_dataframe[['logFC', 'AveExpr', 'P.Value', 'adj.P.Val']]

        # Set index name
        self.limma_dataframe.index.name = 'GeneSymbol'

        # Print
        print 'Done!'

#######################################################
########## 6. Get Differential Expression #############
#######################################################

    def display_deg(self, display_type='table', pvalue_column = 'P.Value', pvalue_threshold = 0.05):

        # Display table
        if display_type == 'table':

            # Print QGrid
            return qgrid.show_grid(self.limma_dataframe)

        elif display_type == 'plot':

            # Filter
            limma_data = {'signif': self.limma_dataframe[self.limma_dataframe[pvalue_column] <= pvalue_threshold],
              'not_signif': self.limma_dataframe[self.limma_dataframe[pvalue_column] > pvalue_threshold]}

            # Create traces
            trace0 = go.Scatter(
                x = limma_data['not_signif']['AveExpr'],
                y = limma_data['not_signif']['logFC'],
                mode = 'markers',
                name = 'Not Significant',
                text = ['<br>'.join([index, 'logFC='+'{:f}'.format(rowData['logFC']), 'p='+'{:.2E}'.format(rowData['P.Value']), 'FDR='+'{:.2E}'.format(rowData['adj.P.Val'])]) for index, rowData in limma_data['not_signif'].iterrows()],
                hoverinfo = 'text+name',
                marker = dict(
                    color = 'black'
                )
                
            )
            trace1 = go.Scatter(
                x = limma_data['signif']['AveExpr'],
                y = limma_data['signif']['logFC'],
                mode = 'markers',
                name = '{pvalue_column} < {pvalue_threshold}'.format(**locals()),
                text = ['<br>'.join([index, 'logFC='+'{:f}'.format(rowData['logFC']), 'p='+'{:.2E}'.format(rowData['P.Value']), 'FDR='+'{:.2E}'.format(rowData['adj.P.Val'])]) for index, rowData in limma_data['signif'].iterrows()],
                hoverinfo = 'text+name',
                    marker = dict(
                    color = 'red'
                )
            )
            layout = go.Layout(
                hovermode = 'closest',
                xaxis = dict(
                    title='Average Expression'
                ),
                yaxis = dict(
                    title='logFoldChange'
                )
            )

            data = [trace0, trace1]
            fig = go.Figure(data=data, layout=layout)
            return(iplot(fig))

        elif display_type == 'plotgl':

            # Create traces
            trace = go.Scattergl(
                x = self.limma_dataframe['AveExpr'],
                y = self.limma_dataframe['logFC'],
                mode = 'markers',
                text = ['<br>'.join([index, 'logFC='+'{:f}'.format(rowData['logFC']), 'p='+'{:.2E}'.format(rowData['P.Value']), 'FDR='+'{:.2E}'.format(rowData['adj.P.Val'])]) for index, rowData in self.limma_dataframe.iterrows()],
                hoverinfo = 'text',
                marker = dict(
                    color = 'black'
                )
            )

            layout = go.Layout(
                hovermode = 'closest',
                xaxis = dict(
                    title='Average Expression'
                ),
                yaxis = dict(
                    title='logFoldChange'
                )
            )

            data = [trace]
            fig = go.Figure(data=data, layout=layout)
            return(iplot(fig))

        else:
            raise ValueError("Parameter display_type must be one in ['table', 'plot', 'plotgl'].")

#######################################################
########## 7. Run Characteristic Direction ############
#######################################################

    def run_characteristic_direction(self, constant_threshold=1e-5):

        # Print
        print 'Running characteristic direction...'

        # Run limma
        self.characteristic_direction_dataframe = pandas2ri.ri2py(r.run_characteristic_direction(pandas2ri.py2ri(self.rawcount_dataframe), pandas2ri.py2ri(self.design_dataframe), constant_threshold=constant_threshold))

        # Set index name
        self.characteristic_direction_dataframe.index.name = 'GeneSymbol'

        # Print
        print 'Done!'

#######################################################
########## 8. Run Enrichr #############################
#######################################################

    def run_enrichr(self, nr_genes=500, gene_set_libraries=['GO_Biological_Process_2015', 'GO_Molecular_Function_2015'], overlapping_genes=False):

        # Print
        print 'Running Enrichr...'

        # Get genesets
        self.genesets = {
            'upregulated': self.characteristic_direction_dataframe.nlargest(nr_genes, 'CD').index.tolist(),
            'downregulated': self.characteristic_direction_dataframe.nsmallest(nr_genes, 'CD').index.tolist(),
        }

        # Add results
        self.enrichr_links = {}
        self.enrichr_results = {}

        # Loop through genesets
        for geneset_label, geneset in self.genesets.iteritems():

            # Add list
            ENRICHR_URL = 'http://amp.pharm.mssm.edu/Enrichr/addList'
            genes_str = '\n'.join(geneset)
            payload = {
                'list': (None, genes_str),
            }
            response = requests.post(ENRICHR_URL, files=payload)
            if not response.ok:
                raise Exception('Error analyzing gene list')
            data = json.loads(response.text)
            user_list_id = data['userListId']
            short_id = data['shortId']

            # Get results
            ENRICHR_URL = 'http://amp.pharm.mssm.edu/Enrichr/enrich'
            query_string = '?userListId=%s&backgroundType=%s'
            result_list = []
            for gene_set_library in gene_set_libraries:
                sleep(0.5)
                response = requests.get(
                    ENRICHR_URL + query_string % (user_list_id, gene_set_library)
                 )
                if not response.ok:
                    raise Exception('Error fetching enrichment results')

                data = json.loads(response.text)
                resultDataframe = pd.DataFrame(data[gene_set_library], columns=['rank', 'term_name', 'pvalue', 'zscore', 'combined_score', 'overlapping_genes', 'FDR', 'old_pvalue', 'old_FDR'])
                selectedColumns = ['term_name','zscore','combined_score','FDR'] if not overlapping_genes else ['term_name','zscore','combined_score','FDR', 'overlapping_genes']
                resultDataframe = resultDataframe.loc[:,selectedColumns]
                resultDataframe['gene_set_library'] = gene_set_library
                result_list.append(resultDataframe)

            # Get results
            mergedResultDataframe = pd.concat(result_list).sort_values('FDR')

            # Add results
            self.enrichr_links[geneset_label] = "http://amp.pharm.mssm.edu/Enrichr/enrich?dataset={short_id}".format(**locals())
            self.enrichr_results[geneset_label] = mergedResultDataframe

        # Print
        print 'Done!'

#######################################################
########## 9. Display Enrichment Results ##############
#######################################################

    def display_enrichment_results(self, display_type='barchart', pvalue_threshold=False, geneset=None):

        # Barchart
        if display_type == 'barchart':

            # Get figures
            figures = []

            # Loop through results
            for geneset_label, enrichment_results in self.enrichr_results.iteritems():
                
                # Get subset and reverse
                enrichment_plot_dataframe = enrichment_results.iloc[:20].iloc[::-1]
                
                # Get logP
                enrichment_plot_dataframe['logFDR'] = -np.log10(enrichment_plot_dataframe['FDR'])
                    
                # Make trace
                trace = go.Bar(
                    x = enrichment_plot_dataframe['logFDR'],
                    y = enrichment_plot_dataframe['term_name'],
                    orientation = 'h',
                    text = [rowData['term_name']+'<br>'+rowData['gene_set_library'].replace('_', ' ').title()+'<br>FDR = '+'{:.2E}'.format(rowData['FDR']) + '<br>Z Score = ' + '{:0.2f}'.format(rowData['zscore'])  + '<br>Combined Score = ' + '{:0.2f}'.format(rowData['combined_score']) for index, rowData in enrichment_plot_dataframe.iterrows()],
                    hoverinfo = 'text',
                    marker=go.Marker(
                        cmax=max(enrichment_plot_dataframe['combined_score']),
                        cmin=0,
                        color=enrichment_plot_dataframe['combined_score'],
                        colorbar=go.ColorBar(
                            title='Combined<br>Score'
                        ),
                        colorscale='Blues' if geneset_label == 'downregulated' else 'Reds'
                    ),
                )

                # Set line
                if pvalue_threshold:
                    shapes=[
                        {
                            'type': 'line',
                            'xref': 'x',
                            'yref': 'paper',
                            'x0': -np.log10(pvalue_threshold),
                            'y0': 0,
                            'x1': -np.log10(pvalue_threshold),
                            'y1': 1,
                            'line': {
                                'color': '#8B0000',
                                'width': 2,
                            },
                        }
                    ]
                else:
                    shapes = []
                
                # Make layout
                layout = go.Layout(
                    margin=go.Margin(
                        l=500,
                        r=50,
                        b=100,
                        t=100,
                        pad=4
                    ),
                    xaxis = dict(
                        title='logFDR'
                    ),
                    yaxis = dict(
                        title='Term'
                    ),
                    title='Enrichment of {geneset_label} genes'.format(**locals()),
                    shapes=shapes
                )
                
                # Make figure
                fig = go.Figure(data=[trace], layout=layout)

                # Plot
                iplot(fig)
                
        # Table
        elif display_type == 'table':

            # Specify geneset
            if not geneset:
                raise ValueError("Please specify 'geneset' parameter: must be in ['upregulated', 'downregulated'].")

            # Print QGrid
            return qgrid.show_grid(self.enrichr_results[geneset].set_index('term_name'))

        # Links
        elif display_type == 'links':

            # Create link dataframe
            link_dataframe = pd.DataFrame.from_dict(self.enrichr_links, orient='index').rename(columns={0:'URL'})

            # Fix URL
            link_dataframe['URL'] = ['<a href="{x}">{x}</a>'.format(**locals()) for x in link_dataframe['URL']]

            # Add index label
            link_dataframe.index.name = 'Geneset'

            # Return HTML
            return HTML(link_dataframe.to_html(escape=False))

#######################################################
########## 10. Run L1000CDS2 ##########################
#######################################################

    def run_l1000cds2(self, mimic=[True, False]):

        # Print
        print 'Running L1000CDS2...'

        # Add results
        # self.characteristic_direction_dataframe = self.rawcount_dataframe.rename(columns={self.rawcount_dataframe.columns[0]:'CD'})
        self.l1000cds2_links = {}
        self.l1000cds2_results = {}

        # Define result dataframe and list
        resultSignatureDataframe = pd.DataFrame()
        linkList = []

        # Loop through aggravate
        for aggravate in mimic:

            # Label
            direction = 'mimic' if aggravate else 'reverse'
            
            # Set data
            data = {"genes": self.characteristic_direction_dataframe.index.tolist(), "vals":self.characteristic_direction_dataframe['CD'].tolist()}
            data['genes'] = [x.upper() for x in data['genes']]
            
            # Set configuration
            config = {"aggravate":aggravate, "searchMethod":"CD", "share":True, "combination":False, "db-version":"latest"}
            payload = {"data":data,"config":config}
            headers = {'content-type':'application/json'}
            
            # Perform request
            r = requests.post('http://amp.pharm.mssm.edu/L1000CDS2/query',data=json.dumps(payload),headers=headers)
            resCD= r.json()
            
            # Add results
            self.l1000cds2_links[direction] = 'http://amp.pharm.mssm.edu/L1000CDS2/#/result/' + resCD['shareId']
            self.l1000cds2_results[direction] = pd.DataFrame(resCD['topMeta']).drop('overlap', axis=1).replace('-666', np.nan)

        # Print
        print 'Done!'

#######################################################
########## 11. Display L1000CDS2 Results ##############
#######################################################

    def display_l1000cds2_results(self, display_type='barchart', direction=None):

        # Group bt
        group_by='pert_desc'

        # Display barchart
        if display_type == 'barchart':

            # Define figures
            figures = []

            # Loop through aggravate
            for direction, l1000cds2_result_dataframe in self.l1000cds2_results.iteritems():

                # Get summary dataframe
                summary_dataframe = l1000cds2_result_dataframe[[group_by,'score']].groupby(group_by).agg(['mean', 'count'])['score'].sort_values('count').replace('', np.nan).dropna()

                # Make trace
                trace = go.Bar(
                    x = summary_dataframe['count'],
                    y = summary_dataframe.index,
                    orientation = 'h',
                    text = 'hello',
                    hoverinfo = 'text',
                    marker=go.Marker(
                        color=summary_dataframe['mean'],
                        colorbar=go.ColorBar(
                            title='Average Score'
                        ),
                        colorscale='Blues' if direction == 'reverse' else 'Reds'
                    ),
                )

                # Make layout
                layout = go.Layout(
                    margin=go.Margin(
                            l=200,
                            r=50,
                            b=100,
                            t=100,
                            pad=4
                        ),
                    xaxis = dict(
                        title='Count (times within top 50)'
                    ),
                    yaxis = dict(
                        title='Perturbagen'
                    ),
                    title = direction
                )

                # Make figure
                fig = go.Figure(data=[trace], layout=layout)

                # Plot
                iplot(fig)

        # Display table
        elif display_type == 'table':

            # Specify geneset
            if not direction:
                raise ValueError("Please specify 'direction' parameter: must be in ['mimic', 'reverse'].")

            # Print QGrid
            return qgrid.show_grid(self.l1000cds2_results[direction].set_index('cell_id'))

        # Display links
        elif display_type == 'links':

            # Create link dataframe
            link_dataframe = pd.DataFrame.from_dict(self.l1000cds2_links, orient='index').rename(columns={0:'URL'})

            # Fix URL
            link_dataframe['URL'] = ['<a href="{x}">{x}</a>'.format(**locals()) for x in link_dataframe['URL']]

            # Add index label
            link_dataframe.index.name = 'Direction'

            # Return HTML
            return HTML(link_dataframe.to_html(escape=False))


#######################################################
########## .  ##################################
#######################################################

#################################################################
#################################################################
########## 2. Fetch Data ########################################
#################################################################
#################################################################

def fetch_dataset(dataset_url, dataset_type):
    pass