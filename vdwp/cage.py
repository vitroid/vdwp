# coding: utf-8
import CageIntegral.physconst as pc
from logging import getLogger
from functools import lru_cache
import CageIntegral.histo as histo
import CageIntegral.histo2f as histo2f


# @lru_cache
def EncagingFE(temperatures, guest, stericterm):
    logger = getLogger()
    f_c = dict()
    for cage in (12, 14, 15, 16):
        histofile = f"data/{guest}.ice{cage}.histo"
        histogram = histo.loadAHisto(open(histofile))

        # f_c ######################
        if histogram is not None:

            f_c[cage] = (histo2f.fvalue(histogram, temperatures) + stericterm)
            logger.debug(f"{cage}: {f_c[cage]}")
    return f_c
