import logging
from _helpers import mock_snakemake, configure_logging

logger = logging.getLogger(__name__)

if __name__ == "__main__":
    if "snakemake" not in globals():
        snakemake = mock_snakemake("log_test")

    configure_logging(snakemake, logger=logger, skip_handlers=False)
    logger.info(logger.handlers)
    print(logger.handlers)
    logger.info("Test log message")
