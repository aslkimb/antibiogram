
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
import datetime
import os

# --- Configuration ---
ISOLATES_FILE = 'Data/Isolates per specimen.csv'
AST_FILE = 'Data/AST.csv'
OUTPUT_FILE = 'Antibiogram_Report.html'

# Gram Stain Dictionary (You might need to expand this based on your data)
GRAM_NEG = [
    'Escherichia coli', 'Klebsiella pneumoniae', 'Pseudomonas aeruginosa', 
    'Proteus mirabilis', 'Klebsiella aerogenes', 'Enterobacter cloacae',
    'Citrobacter freundii', 'Acinetobacter Baumanni', 'Morganella morgani',
    'Salmonella Species', 'Shigella Species', 'Neisseria gonorrhoeae', 
    'Haemophilus influenzae', 'Burkholderia cepacia', 'Stenotrophomonas maltophilia',
    'Serratia marcescens', 'Escherichia coli O157'
]

GRAM_POS = [
    'Staphylococcus aureus', 'Staphylococcus epidermidis', 'Staphylococcus saprophyticus',
    'Staphylococcus hemolyticus', 'Enterococcus faecalis', 'Enterococcus faecium',
    'Streptococcus pneumoniae', 'Streptococcus pyogenes', 'Streptococcus agalactiae',
    'Streptococcus viridians', 'Coagulase negative Staphylococcus', 'Listeria monocytogenes',
    'Staphylococcus lentus', 'Staphylococcus sciuri', 'Aerococcus viridans',
    'Micrococcus species', 'Bacillus Species'
]
    # Fungi are typically excluded from bacterial antibiograms but we can include them in a separate section if needed.
# For this "World Class" view, we will focus on Bacteria for the main heatmaps, or have a separate one for Candida if requested.
# Let's check the data. There are many Candida. Let's make a separate group for Fungi/Yeast.
FUNGI = [
    'Candida albicans', 'Candida glabrata', 'Candida tropicalis', 'Candida parapsilosis', 
    'Candida krusei', 'Candida auris', 'Cryptococcus neoformans', 'Candida guillermondii',
    'Candida dubliniensis', 'Cryptococcus laurentii', 'Candida Lusitaniae', 'Trichosporon asahii'
]

# Known Antifungals (to filter out noise in Fungal section)
ANTIFUNGALS = [
    'Fluconazole', 'Voriconazole', 'Caspofungin', 'Micafungin', 'Flucytosine', 
    'Amphotericin B', 'Itraconazole', 'Posaconazole', 'Anidulafungin'
]


def load_and_clean_data():

    """Loads and cleans the datasets."""
    print("Loading data...")
    try:
        isolates_df = pd.read_csv(ISOLATES_FILE)
        ast_df = pd.read_csv(AST_FILE)
    except FileNotFoundError as e:
        print(f"Error: {e}")
        return None, None

    # --- Preprocessing Isolates ---
    # Parse dates
    isolates_df['Created on'] = pd.to_datetime(isolates_df['Created on'], errors='coerce')
    
    # Rename columns for consistency if needed (based on file view)
    # ISOLATES: Specimen,Organism,AST,Patient,Comment,Relevance,Created on,Sample Type
    # AST: Antimicrobial,MIC,Interpretation,Notes,Isolates,Specimen,Staph CLSI Tier,Drug on Market,Enterobact CLSI Tier
    
    # We need to link AST to Isolates. 
    # The 'Isolates' column in AST.csv looks like '224975030-400' and 'Specimen' column in AST.csv looks like '224975030'.
    # The 'Specimen' column in Isolates.csv looks like '23399301001'.
    # Let's verify the link key. 
    # In AST.csv, 'Specimen' seems to match 'Specimen' in Isolates.csv?
    # Let's assume 'Specimen' is the join key.
    
    # Convert Specimen to string and handle potential float conversion artifacts (e.g. "123.0")
    isolates_df['Specimen'] = isolates_df['Specimen'].astype(str).str.replace(r'\.0$', '', regex=True)
    ast_df['Specimen'] = ast_df['Specimen'].astype(str).str.replace(r'\.0$', '', regex=True)
    
    # --- Deduplication (CLSI M39) ---
    print("Deduplicating isolates (First isolate per patient per species)...")
    
    # Sort by date
    isolates_df = isolates_df.sort_values(by=['Patient', 'Organism', 'Created on'])
    
    # Drop duplicates, keeping the first
    unique_isolates_df = isolates_df.drop_duplicates(subset=['Patient', 'Organism'], keep='first')
    
    print(f"Original Isolates: {len(isolates_df)}")
    print(f"Unique Isolates (Deduplicated): {len(unique_isolates_df)}")
    
    return unique_isolates_df, ast_df

def calculate_susceptibility(unique_isolates_df, ast_df):
    """Calculates % Susceptible for each Organism-Antibiotic pair."""
    print("Calculating susceptibility...")
    
    # Merge AST data with unique isolates
    # We only want AST data for the isolates we kept
    merged_df = pd.merge(ast_df, unique_isolates_df[['Specimen', 'Organism', 'Sample Type']], on='Specimen', how='inner')
    
    if len(merged_df) == 0:
        print("WARNING: Merge resulted in 0 rows! Check Specimen IDs.")
        return pd.DataFrame() # Return empty but valid DF to avoid crash later if we handle it


    
    # Normalize Interpretation
    # Common values: 'Susceptible', 'Resistant', 'Intermediate'
    # CLSI M39: Usually exclude 'Intermediate' or treat as 'Resistant'. World-class standard often excludes them from the numerator (S) but includes in denominator (Total).
    # %S = (Number Susceptible / (Number Susceptible + Number Intermediate + Number Resistant)) * 100
    
    merged_df['Interpretation'] = merged_df['Interpretation'].str.strip().str.title()
    
    # Pivot table: Count of S per Org/Drug
    # We need total tested and count susceptible
    
    results = []
    
    organisms = merged_df['Organism'].unique()
    
    for org in organisms:
        org_df = merged_df[merged_df['Organism'] == org]
        
        # Determine Group
        group = 'Other'
        if org in GRAM_NEG: group = 'Gram-Negative'
        elif org in GRAM_POS: group = 'Gram-Positive'
        elif org in FUNGI: group = 'Fungi'
        
        # Isolate count for this organism (from unique_isolates_df to be accurate about prevalence)
        # Note: The merged_df only has isolates that *had* AST done. Some isolates in unique_isolates_df might not have AST? 
        # But for antibiogram, we care about tested isolates.
        # num_isolates = len(unique_isolates_df[unique_isolates_df['Organism'] == org]) # This is prevalence
        num_isolates_tested = len(org_df['Specimen'].unique()) 
        
        if num_isolates_tested < 1: continue

        antibiotics = org_df['Antimicrobial'].unique()
        
        for abx in antibiotics:
            abx_df = org_df[org_df['Antimicrobial'] == abx]
            
            total_tested = len(abx_df)
            if total_tested == 0: continue
            
            susceptible_count = len(abx_df[abx_df['Interpretation'] == 'Susceptible'])
            percent_s = (susceptible_count / total_tested) * 100
            
            results.append({
                'Organism': org,
                'Antibiotic': abx,
                'Group': group,
                'Isolates_Tested': total_tested, 
                'Percent_S': percent_s,
                'Total_Isolates_Of_Org': num_isolates_tested
            })
            
    return pd.DataFrame(results)

def create_heatmap(df, group_name, colorscale='RdYlGn'):
    """Creates a Plotly heatmap for a specific group."""
    if df.empty:
        return None
        
    # Pivot for Heatmap
    # Rows: Organism (with count), Cols: Antibiotic
    
    # Add count to Organism name for display
    df['Org_Label'] = df['Organism'] + " (n=" + df['Total_Isolates_Of_Org'].astype(str) + ")"
    
    heatmap_data = df.pivot(index='Org_Label', columns='Antibiotic', values='Percent_S')
    
    # Fill NaN with None so they show up as empty/grey
    # But Plotly handles NaNs well in heatmaps (transparent).
    
    # Hover Text
    hover_text = []
    for index, row in heatmap_data.iterrows():
        hover_row = []
        for col in heatmap_data.columns:
            val = row[col]
            if pd.isna(val):
                hover_row.append("Not Tested")
            else:
                # Find the N for this specific cell
                # We need to look back at the original df
                # This is a bit slow doing it this way, let's optimize if dataset is huge, but it's small.
                record = df[(df['Org_Label'] == index) & (df['Antibiotic'] == col)]
                if not record.empty:
                    n_tested = record.iloc[0]['Isolates_Tested']
                    hover_row.append(f"{val:.1f}% Susceptible<br>Tested: {n_tested}")
                else:
                    hover_row.append("Error")
        hover_text.append(hover_row)

    fig = go.Figure(data=go.Heatmap(
        z=heatmap_data.values,
        x=heatmap_data.columns,
        y=heatmap_data.index,
        text=hover_text,
        hoverinfo='text',
        colorscale=colorscale,
        zmin=0, zmax=100,
        xgap=1, ygap=1
    ))
    
    fig.update_layout(
        title=f'{group_name} Bacteria Susceptibility (% Susceptible)',
        xaxis_title="Antibiotic",
        yaxis_title="Organism",
        height=max(400, len(heatmap_data)*30), # Adjust height based on rows
        yaxis={'autorange': 'reversed'} # Common for heatmaps to have A at top
    )
    
    return fig

def generate_dashboard():
    unique_isolates, ast_df = load_and_clean_data()
    if unique_isolates is None: return

    susceptibility_df = calculate_susceptibility(unique_isolates, ast_df)
    
    # Filter for significant organisms (optional, but requested "World Class" often filters <30)
    # For this demo, let's keep everything but maybe sort?
    # Or separate Low Count isolates.
    # Let's keep all for visibility but add a warning in the HTML if N < 30.
    
    # Split by Group
    gram_neg_df = susceptibility_df[susceptibility_df['Group'] == 'Gram-Negative']
    gram_pos_df = susceptibility_df[susceptibility_df['Group'] == 'Gram-Positive']
    
    # Fungi Group - Filter out non-antifungals
    # Some datasets might map antibacterials to Candida erroneously.
    fungi_df = susceptibility_df[susceptibility_df['Group'] == 'Fungi']
    if not fungi_df.empty:
        fungi_df = fungi_df[fungi_df['Antibiotic'].isin(ANTIFUNGALS)]
    
    # Create Figures for General Overview
    fig_gn = create_heatmap(gram_neg_df, "Gram-Negative")
    fig_gp = create_heatmap(gram_pos_df, "Gram-Positive")
    fig_fungi = create_heatmap(fungi_df, "Fungal")

    # --- Specimen Type Analysis ---
    # Find Top 5 Specimen Types by isolate count (to ensure Blood is included if present, often rank #4)
    top_specimens = unique_isolates['Sample Type'].value_counts().head(5).index.tolist()
    
    specimen_sections_html = ""
    
    for specimen in top_specimens:
        specimen_clean_name = specimen.replace(" ", "_")
        
        # Filter data for this specimen
        # We need to filter the susceptibility_df? No, we need to recalculate susceptibility for this subset.
        # But calculate_susceptibility takes (unique_isolates_df, ast_df).
        # So we just filter unique_isolates first!
        
        specimen_isolates = unique_isolates[unique_isolates['Sample Type'] == specimen]
        if len(specimen_isolates) < 1: continue
        
        specimen_susceptibility = calculate_susceptibility(specimen_isolates, ast_df)
        
        if specimen_susceptibility.empty: continue
        
        # Determine main group for this specimen (usually we care about all, but splitting heatmaps is better)
        # For simplicity in this section, let's create one Combined heatmap for top pathogens in this specimen
        # or split Gram Neg/Pos again. Splitting is cleaner.
        
        spec_gn = specimen_susceptibility[specimen_susceptibility['Group'] == 'Gram-Negative']
        spec_gp = specimen_susceptibility[specimen_susceptibility['Group'] == 'Gram-Positive']
        
        fig_spec_gn = create_heatmap(spec_gn, f"{specimen} - Gram-Negative")
        fig_spec_gp = create_heatmap(spec_gp, f"{specimen} - Gram-Positive")
        
        specimen_html = f"""
        <div class="card">
            <div class="card-header bg-light">
                Specimen Analysis: {specimen} <span class="badge bg-secondary">{len(specimen_isolates)} Isolates</span>
            </div>
            <div class="card-body">
                {fig_spec_gn.to_html(full_html=False, include_plotlyjs=False) if fig_spec_gn else '<p class="text-muted">No Gram-Negative Data</p>'}
                <hr>
                {fig_spec_gp.to_html(full_html=False, include_plotlyjs=False) if fig_spec_gp else '<p class="text-muted">No Gram-Positive Data</p>'}
            </div>
        </div>
        """
        specimen_sections_html += specimen_html


    # --- Encode Logo ---
    logo_path = r"C:/Users/LabAdmin/.gemini/antigravity/brain/2b5c8207-f539-4bdb-a291-648596afb150/uploaded_media_1769448383060.jpg"
    logo_html = ""
    try:
        import base64
        with open(logo_path, "rb") as image_file:
            encoded_string = base64.b64encode(image_file.read()).decode('utf-8')
            logo_html = f'<img src="data:image/jpeg;base64,{encoded_string}" alt="Nile International Hospital Logo" style="max-height: 120px; margin-bottom: 1rem;">'
    except Exception as e:
        print(f"Warning: Could not load logo: {e}")

    # --- Generate HTML ---
    # We will write raw HTML and embed the Plotly divs.
    
    # Branding Colors extracted from Logo
    COLOR_NAVY = "#1F3A60" # Dark Blue text
    COLOR_CYAN = "#26C6DA" # The 'I' figure and wave
    
    html_content = f"""
    <!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>Antibiogram - Nile International Hospital</title>
        <!-- Bootstrap 5 -->
        <link href="https://cdn.jsdelivr.net/npm/bootstrap@5.3.0/dist/css/bootstrap.min.css" rel="stylesheet">
        <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
        <style>
            body {{ background-color: #f4f6f9; font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif; color: #333; }}
            
            /* Header: Clean White to display Logo */
            .header {{ background-color: white; padding: 2rem 0; margin-bottom: 2rem; border-bottom: 5px solid {COLOR_CYAN}; box-shadow: 0 2px 10px rgba(0,0,0,0.05); }}
            
            /* Typography */
            h1, h2, h3 {{ color: {COLOR_NAVY}; }}
            .brand-subtitle {{ color: {COLOR_CYAN}; font-weight: 500; letter-spacing: 1px; text-transform: uppercase; font-size: 1.1rem; }}
            
            /* Cards */
            .card {{ border: none; box-shadow: 0 4px 6px rgba(0,0,0,0.05); margin-bottom: 2rem; border-radius: 8px; overflow: hidden; }}
            .card-header {{ background-color: {COLOR_NAVY}; color: white; border-bottom: none; font-weight: 600; font-size: 1.25rem; padding: 1rem 1.5rem; }}
            .card-header .badge {{ background-color: {COLOR_CYAN} !important; color: {COLOR_NAVY}; }}
            
            /* Stats */
            .stat-box {{ text-align: center; padding: 1.5rem; background: white; border-radius: 8px; border: 1px solid #eee; transition: transform 0.2s; }}
            .stat-box:hover {{ transform: translateY(-5px); box-shadow: 0 10px 20px rgba(0,0,0,0.1); }}
            .stat-value {{ font-size: 2.5rem; font-weight: bold; color: {COLOR_NAVY}; }}
            .stat-label {{ color: {COLOR_CYAN}; font-size: 0.85rem; text-transform: uppercase; letter-spacing: 1px; font-weight: 700; }}
            
            /* Footer */
            .footer {{ text-align: center; padding: 3rem 0; color: #6c757d; font-size: 0.9rem; border-top: 1px solid #eee; background-color: white; margin-top: 3rem; }}
            .footer .author {{ font-weight: bold; color: {COLOR_NAVY}; }}
        </style>
    </head>
    <body>

    <div class="header">
        <div class="container text-center">
            {logo_html}
            <p class="brand-subtitle mb-2">– Source of Health –</p>
            <h2 class="fw-light mb-0" style="color: #666;">Cumulative Antibiogram Report</h2>
            <small class="text-muted">Generated: {datetime.datetime.now().strftime("%Y-%m-%d %H:%M")}</small>
        </div>
    </div>

    <div class="container">
        
        <!-- Summary Stats -->
        <div class="row mb-5">
            <div class="col-md-3">
                <div class="card stat-box">
                    <div class="stat-value">{len(unique_isolates)}</div>
                    <div class="stat-label">Unique Isolates</div>
                </div>
            </div>
             <div class="col-md-3">
                <div class="card stat-box">
                    <div class="stat-value">{len(susceptibility_df['Organism'].unique())}</div>
                    <div class="stat-label">Organisms Identified</div>
                </div>
            </div>
             <div class="col-md-3">
                <div class="card stat-box">
                    <div class="stat-value">{int(susceptibility_df['Percent_S'].mean())}%</div>
                    <div class="stat-label">Avg. Susceptibility</div>
                </div>
            </div>
             <div class="col-md-3">
                <div class="card stat-box">
                    <div class="stat-value">{len(susceptibility_df['Antibiotic'].unique())}</div>
                    <div class="stat-label">Antibiotics Tested</div>
                </div>
            </div>
        </div>

        <div class="alert alert-light border d-flex align-items-center mb-5" role="alert" style="border-left: 5px solid {COLOR_NAVY} !important;">
            <div style="color: {COLOR_NAVY}; margin-right: 15px;">
                <svg xmlns="http://www.w3.org/2000/svg" width="32" height="32" fill="currentColor" class="bi bi-info-circle-fill" viewBox="0 0 16 16">
                  <path d="M8 16A8 8 0 1 0 8 0a8 8 0 0 0 0 16zm.93-9.412-1 4.705c-.07.34.029.533.304.533.194 0 .487-.07.686-.246l-.088.416c-.287.346-.92.598-1.465.598-.703 0-1.002-.422-.808-1.319l.738-3.468c.064-.293.006-.399-.287-.47l-.451-.081.082-.381 2.29-.287zM8 5.5a1 1 0 1 1 0-2 1 1 0 0 1 0 2z"/>
                </svg>
            </div>
            <div>
                <strong style="color: {COLOR_NAVY};">Methodology:</strong><br>
                Data follows CLSI M39 guidelines. Only the first isolate per patient per species is included.<br>
                Values represent the percentage of isolates susceptible to the antibiotic.
            </div>
        </div>

        <!-- Gram Negative Section -->
        <div class="card">
            <div class="card-header">
                Gram-Negative Bacteria
            </div>
            <div class="card-body">
                {fig_gn.to_html(full_html=False, include_plotlyjs='cdn') if fig_gn else '<p class="text-center text-muted">No Gram-Negative Data</p>'}
            </div>
        </div>

        <!-- Gram Positive Section -->
        <div class="card">
            <div class="card-header">
                Gram-Positive Bacteria
            </div>
            <div class="card-body">
                {fig_gp.to_html(full_html=False, include_plotlyjs=False) if fig_gp else '<p class="text-center text-muted">No Gram-Positive Data</p>'}
            </div>
        </div>
        
         <!-- Fungi Section -->
        <div class="card">
            <div class="card-header">
                Fungi & Yeasts (Antifungals Only)
            </div>
            <div class="card-body">
                {fig_fungi.to_html(full_html=False, include_plotlyjs=False) if fig_fungi else '<p class="text-center text-muted">No Fungal Data</p>'}
            </div>
        </div>
        
        <h2 class="mt-5 mb-3">Stratification by Specimen Type</h2>
        <p class="text-muted mb-4">Detailed analysis for the top 5 most common specimen types.</p>
        
        {specimen_sections_html}

    </div>

    <div class="footer">
        <p>Prepared by <span class="author">Aslam Kimbugwe</span> - Medical Laboratory Scientist</p>
        <p>&copy; {datetime.datetime.now().year} Nile International Hospital. All Rights Reserved.</p>
    </div>

    </body>
    </html>
    """

    # --- Generate PDF ---
    def generate_pdf():
        print("Generating PDF Report...")
        from fpdf import FPDF
        import os
        
        pdf_file = "Antibiogram_Report.pdf"
        
        class PDF(FPDF):
            def header(self):
                # We will add a custom header on each page if needed, but for now cover page is unique
                pass
                
            def footer(self):
                self.set_y(-15)
                self.set_font('Helvetica', 'I', 8)
                self.set_text_color(128)
                self.cell(0, 10, f'Page {self.page_no()} | Prepared by Aslam Kimbugwe - Nile International Hospital', align='C')

        pdf = PDF()
        pdf.set_auto_page_break(auto=True, margin=15)
        
        # --- Cover Page ---
        pdf.add_page()
        
        # Logo
        try:
            # Save base64 logo to temp file for FPDF or use original path
            # Using original path is safer
            pdf.image(logo_path, x=65, y=30, w=80) 
        except:
            pass
            
        pdf.ln(100)
        
        # Title
        pdf.set_font("Helvetica", "B", 24)
        pdf.set_text_color(31, 58, 96) # Navy #1F3A60
        pdf.cell(0, 10, "NILE INTERNATIONAL HOSPITAL", align='C', new_x="LMARGIN", new_y="NEXT")
        
        pdf.set_font("Helvetica", "I", 14)
        pdf.set_text_color(38, 198, 218) # Cyan #26C6DA
        pdf.cell(0, 10, "- Source of Health -", align='C', new_x="LMARGIN", new_y="NEXT")
        
        pdf.ln(20)
        pdf.set_font("Helvetica", "", 18)
        pdf.set_text_color(50)
        pdf.cell(0, 10, "Cumulative Antibiogram Report", align='C', new_x="LMARGIN", new_y="NEXT")
        
        pdf.set_font("Helvetica", "", 12)
        pdf.cell(0, 10, f"Generated: {datetime.datetime.now().strftime('%Y-%m-%d')}", align='C', new_x="LMARGIN", new_y="NEXT")

        # --- Summary Stats Page ---
        pdf.add_page()
        pdf.set_font("Helvetica", "B", 16)
        pdf.set_text_color(31, 58, 96)
        pdf.cell(0, 10, "Executive Summary", new_x="LMARGIN", new_y="NEXT")
        pdf.ln(5)
        
        pdf.set_font("Helvetica", "", 12)
        pdf.set_text_color(0)
        
        # Simple list of stats
        stats = [
            f"Unique Isolates Processed: {len(unique_isolates)}",
            f"Distinct Organisms Identified: {len(susceptibility_df['Organism'].unique())}",
            f"Average Susceptibility Rate: {int(susceptibility_df['Percent_S'].mean())}%",
            f"Antibiotics Tested: {len(susceptibility_df['Antibiotic'].unique())}"
        ]
        
        for stat in stats:
            pdf.cell(0, 10, chr(149) + " " + stat, new_x="LMARGIN", new_y="NEXT")
            
        pdf.ln(10)
        pdf.set_font("Helvetica", "I", 10)
        pdf.set_text_color(100)
        pdf.multi_cell(0, 6, "Methodology: Data follows CLSI M39 guidelines. Only the first isolate per patient per species is included. Analysis excludes intermediate results from % Susceptible calculation unless otherwise noted.")

        # --- Heatmaps ---
        # Helper to save and add plot
        def add_plot_to_pdf(fig, title):
            if not fig: return
            pdf.add_page()
            pdf.set_font("Helvetica", "B", 16)
            pdf.set_text_color(31, 58, 96)
            pdf.cell(0, 10, title, new_x="LMARGIN", new_y="NEXT")
            pdf.ln(5)
            
            # Save temp image
            temp_img = f"temp_{title.replace(' ', '_')}.png"
            try:
                # Kaleido export
                # Increase scale for quality
                fig.write_image(temp_img, width=1000, height=min(1200, max(600, len(fig.data[0].y)*30)), scale=2)
                
                # Add to PDF
                # Adjust width to fit page (A4 width ~210mm)
                pdf.image(temp_img, x=10, w=190)
                
                os.remove(temp_img)
            except Exception as e:
                pdf.set_font("Courier", "", 10)
                pdf.set_text_color(255, 0, 0)
                pdf.cell(0, 10, f"Error rendering chart: {e}", new_x="LMARGIN", new_y="NEXT")
                print(f"Error exporting image: {e}")

        add_plot_to_pdf(fig_gn, "Gram-Negative Bacteria")
        add_plot_to_pdf(fig_gp, "Gram-Positive Bacteria")
        add_plot_to_pdf(fig_fungi, "Fungi and Yeasts")

        # Top Specimens
        # We need to recreate these figures because we didn't store them in variables outside the loop
        # Optimized: Refactor the loop above to store figures or just re-run loop here?
        # Re-running loop is safer to avoid extensive refactoring of previous code block
        
        pdf.add_page()
        pdf.set_font("Helvetica", "B", 18)
        pdf.cell(0, 10, "Specimen Specific Analysis", align='C', new_x="LMARGIN", new_y="NEXT")
        pdf.ln(10)
        
        top_specimens = unique_isolates['Sample Type'].value_counts().head(5).index.tolist()
        
        for specimen in top_specimens:
             specimen_isolates = unique_isolates[unique_isolates['Sample Type'] == specimen]
             if len(specimen_isolates) < 1: continue
             specimen_susceptibility = calculate_susceptibility(specimen_isolates, ast_df)
             if specimen_susceptibility.empty: continue
             
             # Re-create figures
             spec_gn = specimen_susceptibility[specimen_susceptibility['Group'] == 'Gram-Negative']
             spec_gp = specimen_susceptibility[specimen_susceptibility['Group'] == 'Gram-Positive']
             
             fig_s_gn = create_heatmap(spec_gn, f"{specimen} - Gram-Negative")
             fig_s_gp = create_heatmap(spec_gp, f"{specimen} - Gram-Positive")
             
             if fig_s_gn: add_plot_to_pdf(fig_s_gn, f"{specimen}: Gram-Negative")
             if fig_s_gp: add_plot_to_pdf(fig_s_gp, f"{specimen}: Gram-Positive")

        try:
            pdf.output(pdf_file)
            print(f"PDF Report generated successfully: {os.path.abspath(pdf_file)}")
        except Exception as e:
            print(f"Error saving PDF: {e}")

    # Generate the PDF
    try:
        generate_pdf()
    except Exception as e:
        print(f"Could not generate PDF: {e}")
        import traceback
        traceback.print_exc()

    with open(OUTPUT_FILE, 'w', encoding='utf-8') as f:
        f.write(html_content)

    print(f"Dashboard generated successfully: {os.path.abspath(OUTPUT_FILE)}")

if __name__ == "__main__":
    generate_dashboard()
