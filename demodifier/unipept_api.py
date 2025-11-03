# unipept_api.py
# ---------------------------------------------------------
# Handles communication with the Unipept API.
# Each thread uses its own requests.Session for stability.
# ---------------------------------------------------------

import logging
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

# Configure module logger (level/handlers managed by settings.setup_logging)
logger = logging.getLogger(__name__)

# Retry failed HTTP requests up to 3 times with exponential backoff
# Shared retry configuration for all sessions
retries = Retry(
    total=3,
    backoff_factor=0.5,
    status_forcelist=[429, 500, 502, 503, 504],
    raise_on_status=False,
)

def make_session():
    """
    Create a fresh requests.Session with the standard retry policy.

    Each worker thread should use its own session to avoid cross-thread
    reuse of network connections.

    Returns
    -------
    requests.Session
        A new session object configured with retry logic.
    """
    s = requests.Session()
    s.mount("https://", HTTPAdapter(max_retries=retries))
    return s


def process_peptides(peptides, session):
    """
    Query the Unipept pept2lca endpoint for a list of peptides.

    Parameters
    ----------
    peptides : list of str
        Peptide sequences to query.
    session : requests.Session
        Session object used for API requests.

    Returns
    -------
    list of str
        LCAs (lowest common ancestors) aligned to the input peptide order.
        If the request fails, returns a list of "no response" entries.
    """
    url = "https://api.unipept.ugent.be/api/v1/pept2lca"
    params = {"input[]": peptides, "equate_il": "true"}

    logger.debug(f"Querying Unipept API with peptides: {peptides}")

    try:
        response = session.get(url, params=params, timeout=10)
        response.raise_for_status()
        results = response.json()

        logger.debug(f"API response: {results}")

        lca_map = {r["peptide"]: r.get("taxon_name", "no match") for r in results}
        return [lca_map.get(p, "no match") for p in peptides]

    except Exception as e:
        logger.error(f"API request failed: {e}")
        return ["no response"] * len(peptides)


def get_lcas_for_permutations(permutations, session):
    """
    Look up LCAs for a list of peptide variants in chunks of 100.

    Parameters
    ----------
    permutations : list of str
        Peptide variants to query in batches.
    session : requests.Session
        Session object used for API requests.

    Returns
    -------
    list of str
        LCAs corresponding to each input permutation, maintaining order.
        Returns "no match" for any failed lookups.
    """
    chunk_size = 100
    lca_map = {}

    for i in range(0, len(permutations), chunk_size):
        chunk = permutations[i : i + chunk_size]
        chunk_lcas = process_peptides(chunk, session=session)
        for pep, lca in zip(chunk, chunk_lcas):
            lca_map[pep] = lca

    return [lca_map.get(p, "no match") for p in permutations]
