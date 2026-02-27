import pandas as pd
import statsmodels.api as sm
from statsmodels.formula.api import ols

def run_two_way_anova(adata, gene_name, factor1='Diet', factor2='Genotype'):
    """
    Runs Two-Way ANOVA for a specific gene across two metadata conditions.
    
    Args:
        adata: AnnData object containing the cells of interest.
        gene_name: The gene to test (string).
        factor1: Column name in adata.obs for the first condition (e.g., 'Diet').
        factor2: Column name in adata.obs for the second condition (e.g., 'Genotype').
    """
    # 1. Extract expression data for the gene
    try:
        # Check if gene exists in raw or processed data
        if gene_name in adata.raw.var_names:
            expr_data = adata.raw[:, gene_name].X.toarray().flatten()
        else:
            expr_data = adata[:, gene_name].X.toarray().flatten()
    except KeyError:
        print(f"Gene {gene_name} not found.")
        return None

    # 2. Create a DataFrame for statsmodels
    df = pd.DataFrame({
        'Expression': expr_data,
        factor1: adata.obs[factor1].values,
        factor2: adata.obs[factor2].values
    })

    # 3. Fit the model (Expression ~ Factor1 + Factor2 + Interaction)
    formula = f'Expression ~ C({factor1}) + C({factor2}) + C({factor1}):C({factor2})'
    model = ols(formula, data=df).fit()
    
    # 4. Perform ANOVA
    anova_table = sm.stats.anova_lm(model, typ=2)
    
    return anova_table

# Example Usage:
# result = run_two_way_anova(goblet_adata, 'Muc2', factor1='Diet', factor2='Condition')
# print(result)