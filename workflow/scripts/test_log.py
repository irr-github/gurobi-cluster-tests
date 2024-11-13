import sys  # , os
import logging, traceback
from _helpers import mock_snakemake, configure_logging

logger = logging.getLogger(__name__)

if "snakemake" not in globals():
    snakemake = mock_snakemake("log_test")
    configure_logging(snakemake, logger=logger, skip_handlers=False)
    logger.info("Test log message")
