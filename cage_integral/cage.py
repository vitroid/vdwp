# coding: utf-8
import vdwp.physconst as pc
from logging import getLogger
from functools import lru_cache
import vdwp.histo as histo
import vdwp.histo2f as histo2f


# この関数は、温度とゲストの種類を受け取り、ゲストのエンカージングポテンシャルを計算します。
# 温度はリストで指定され、ゲストは文字列で指定されます。
# sterictermはゲストの立体反発項です。


# @lru_cache
def EncagingFE(temperatures, guest, stericterm):
    logger = getLogger()
    f_c = dict()
    for cage in (12, 14, 15, 16):
        histofile = f"data/{guest}.ice{cage}.histo"
        histogram = histo.loadAHisto(open(histofile))

        # f_c ######################
        if histogram is not None:

            f_c[cage] = histo2f.fvalue(histogram, temperatures) + stericterm
            logger.debug(f"{cage}: {f_c[cage]}")
    return f_c
