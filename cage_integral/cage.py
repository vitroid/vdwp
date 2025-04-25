# coding: utf-8
from logging import getLogger
import loader
import histo2f


# この関数は、温度とゲストの種類を受け取り、ゲストのエンカージングポテンシャルを計算します。
# 温度はリストで指定され、ゲストは文字列で指定されます。
# sterictermはゲストの立体反発項です。


def EncagingFE(temperatures, guest, stericterm):
    logger = getLogger()
    f_c = dict()
    for cage in (12, 14, 15, 16):
        histofile = f"{guest}.{cage}hedra.histo"
        histogram = loader.loadAHisto(open(histofile))

        # f_c ######################
        if histogram is not None:

            f_c[cage] = histo2f.fvalue(histogram, temperatures) + stericterm
            logger.debug(f"{cage}: {f_c[cage]}")
    return f_c
