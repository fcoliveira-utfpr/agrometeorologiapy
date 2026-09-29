"""Classificação climática: Köppen-Geiger, Thornthwaite (1948), Camargo (1991)
mod. Maluf (2000) e zonas de vida de Holdridge.

Cada método tem duas versões:

- `classificacao_<metodo>(...)`: um local (listas/arrays de 12 meses).
- `classificacao_<metodo>_grade(...)`: muitos locais de uma vez (ex.: todos os
  pixels de um raster), com os 12 meses no primeiro eixo, formato (12, ...).

Locais sem dado (NaN em qualquer mês) recebem id 0 e classe ''.
"""

import numpy as np

from .balanco_hidrico import balanco_hidrico_climatologico_grade

__all__ = [
    "classificacao_koppen",
    "classificacao_koppen_grade",
    "classificacao_thornthwaite",
    "classificacao_thornthwaite_grade",
    "classificacao_camargo",
    "classificacao_camargo_grade",
    "classificacao_holdridge",
    "classificacao_holdridge_grade",
]


# ---------------------------------------------------------------------------
# Utilitários
# ---------------------------------------------------------------------------
def _mensal(x, nome):
    x = np.asarray(x, dtype=float)
    if x.shape[0] != 12:
        raise ValueError(f"{nome} deve ter 12 meses no primeiro eixo (formato (12, ...)).")
    return x


def _espacial(x, forma):
    return np.broadcast_to(np.asarray(x, dtype=float), forma)


def _faixas(x, limites, codigos):
    """Código da faixa [limites[i], limites[i+1]) de x; NaN -> 0."""
    idx = np.digitize(np.nan_to_num(x, nan=0.0), limites)
    out = np.asarray(codigos, dtype=np.uint8)[idx]
    return np.where(np.isfinite(x), out, 0).astype(np.uint8)


def _rotulos(codigo, legenda):
    """Array de códigos inteiros -> array de strings (0 -> '')."""
    tabela = np.array([''] + [legenda[k] for k in range(1, max(legenda) + 1)], dtype=object)
    return np.asarray(tabela[codigo], dtype=object)


def _escalar(saida):
    """Converte a saída de uma função _grade de um único local em escalares."""
    return {k: (v.item() if isinstance(v, np.ndarray) else v) for k, v in saida.items()}


def _um_local(*mensais):
    for x in mensais:
        if x is not None and np.ndim(x) != 1:
            raise ValueError("Para vários locais use a versão _grade.")


# Estações astronômicas ponderadas (hemisfério sul; índice 0 = janeiro):
# verão = 1/3 dez + jan + fev + 2/3 mar; outono = 1/3 mar + abr + mai + 2/3 jun;
# inverno = 1/3 jun + jul + ago + 2/3 set; primavera = 1/3 set + out + nov + 2/3 dez.
# No hemisfério norte, verão <-> inverno e outono <-> primavera.
_PESOS_SUL = {
    'verao': {11: 1 / 3, 0: 1, 1: 1, 2: 2 / 3},
    'outono': {2: 1 / 3, 3: 1, 4: 1, 5: 2 / 3},
    'inverno': {5: 1 / 3, 6: 1, 7: 1, 8: 2 / 3},
    'primavera': {8: 1 / 3, 9: 1, 10: 1, 11: 2 / 3},
}
_OPOSTA = {'verao': 'inverno', 'inverno': 'verao', 'outono': 'primavera', 'primavera': 'outono'}


def _pesos(estacao):
    w = np.zeros(12)
    for m, p in _PESOS_SUL[estacao].items():
        w[m] = p
    return w


def _soma_estacao(x, estacao, norte):
    """Soma ponderada de x (12, ...) na estação, trocando de hemisfério onde norte."""
    extra = (1,) * (x.ndim - 1)
    sul = (x * _pesos(estacao).reshape((12,) + extra)).sum(0)
    nor = (x * _pesos(_OPOSTA[estacao]).reshape((12,) + extra)).sum(0)
    return np.where(norte, nor, sul)


# ---------------------------------------------------------------------------
# Köppen-Geiger (Alvares et al., 2013; sazonalidade de Kottek et al., 2006)
# ---------------------------------------------------------------------------
KOPPEN_CLASSES = {
    1: "Af", 2: "Am", 3: "As", 4: "Aw",
    5: "BSh", 6: "BSk", 7: "BWh", 8: "BWk",
    9: "Cfa", 10: "Cfb", 11: "Cfc",
    12: "Csa", 13: "Csb", 14: "Csc",
    15: "Cwa", 16: "Cwb", 17: "Cwc",
    18: "Dfa", 19: "Dfb", 20: "Dfc", 21: "Dfd",
    22: "Dsa", 23: "Dsb", 24: "Dsc", 25: "Dsd",
    26: "Dwa", 27: "Dwb", 28: "Dwc", 29: "Dwd",
    30: "ET", 31: "EF",
}

# Semestres out-mar e abr-set (índice 0 = janeiro).
_OUT_MAR = [9, 10, 11, 0, 1, 2]
_ABR_SET = [3, 4, 5, 6, 7, 8]


def classificacao_koppen_grade(T, P, lat):
    """
    Classificação climática de Köppen-Geiger para muitos locais, pela
    metodologia de Alvares et al. (2013), com a sazonalidade f/s/w dos
    grupos C e D de Kottek et al. (2006).

    O verão é o semestre out-mar no hemisfério sul e abr-set no norte.

    Parâmetros
    ----------
    T : array_like, formato (12, ...)
        Temperatura média mensal do ar, em °C (janeiro a dezembro).
    P : array_like, formato (12, ...)
        Precipitação mensal, em mm/mês.
    lat : float ou array_like
        Latitude, em graus (negativa no hemisfério sul); escalar ou com o
        formato das dimensões espaciais.

    Retorna
    -------
    dict de numpy.ndarray
        'id' (1-31, ver KOPPEN_CLASSES; 0 = sem dado) e 'classe' (ex.: 'Cfa').
    """
    T = _mensal(T, "T")
    P = _mensal(P, "P")
    if T.shape != P.shape:
        raise ValueError("T e P devem ter o mesmo formato.")
    forma = T.shape[1:]
    norte = _espacial(lat, forma) > 0
    valido = np.all(np.isfinite(T), axis=0) & np.all(np.isfinite(P), axis=0)

    tann = T.mean(0)
    tcold = T.min(0)
    thot = T.max(0)
    tmon10 = (T > 10).sum(0)
    rann = P.sum(0)
    rdry = P.min(0)

    om, aset = P[_OUT_MAR], P[_ABR_SET]
    p_verao = np.where(norte, aset.sum(0), om.sum(0))
    p_inverno = np.where(norte, om.sum(0), aset.sum(0))
    psdry = np.where(norte, aset.min(0), om.min(0))
    pwdry = np.where(norte, om.min(0), aset.min(0))
    pswet = np.where(norte, aset.max(0), om.max(0))
    pwwet = np.where(norte, om.max(0), aset.max(0))

    with np.errstate(divide="ignore", invalid="ignore"):
        pct_verao = np.where(rann > 0, p_verao / rann, 0.0)
        pct_inverno = np.where(rann > 0, p_inverno / rann, 0.0)

    # Limiar de aridez (grupo B), em cm: 2T + 14, 2T (chuva de inverno) ou 2T + 28 (de verão).
    p_lim = 2 * tann + 14
    p_lim = np.where(pct_inverno >= 0.7, 2 * tann, p_lim)
    p_lim = np.where(pct_verao >= 0.7, 2 * tann + 28, p_lim)

    grupo_a = tcold >= 18
    grupo_c = (thot > 10) & (tcold > -3) & (tcold < 18)
    grupo_d = (thot > 10) & (tcold <= -3)
    grupo_e = thot <= 10
    grupo_b = rann < 10 * p_lim

    # A: f, m, s, w
    am_lim = 100 - rann / 25
    a_seco = grupo_a & (rdry < 60) & (rdry < am_lim)

    # C e D: s (verão seco), w (inverno seco) e f (nem s nem w)
    is_s = (psdry < pwdry) & (pwwet > 3 * psdry) & (psdry < 40)
    is_w = (pwdry < psdry) & (pswet > 10 * pwdry)
    is_f = ~is_s & ~is_w
    quente = thot >= 22
    temperado = ~quente & (tmon10 >= 4)
    curto = ~quente & ~temperado & (tmon10 >= 1)
    muito_frio = tcold < -38

    cod = np.zeros(forma, dtype=np.uint8)
    regras = [
        (grupo_a & (rdry >= 60), 1),
        (grupo_a & (rdry < 60) & (rdry >= am_lim), 2),
        (a_seco & (psdry < pwdry), 3),
        (a_seco & (psdry >= pwdry), 4),
    ]
    for grupo, base in ((grupo_c, 9), (grupo_d, 18)):
        passo = 3 if base == 9 else 4
        for k, saz in enumerate((is_f, is_s, is_w)):
            b = base + k * passo
            regras += [(grupo & saz & quente, b), (grupo & saz & temperado, b + 1),
                       (grupo & saz & curto, b + 2)]
            if base == 18:
                regras.append((grupo & saz & muito_frio, b + 3))
    regras += [
        (grupo_e & (thot > 0), 30),
        (grupo_e & (thot <= 0), 31),
        # B por último: sobrepõe qualquer outro grupo.
        (grupo_b & (rann >= 5 * p_lim) & (tann >= 18), 5),
        (grupo_b & (rann >= 5 * p_lim) & (tann < 18), 6),
        (grupo_b & (rann < 5 * p_lim) & (tann >= 18), 7),
        (grupo_b & (rann < 5 * p_lim) & (tann < 18), 8),
    ]
    for cond, c in regras:
        cod[cond] = c
    cod[~valido] = 0
    return {'id': cod, 'classe': _rotulos(cod, KOPPEN_CLASSES).astype(str)}


def classificacao_koppen(T, P, lat):
    """
    Classificação climática de Köppen-Geiger de um local (Alvares et al.,
    2013; Kottek et al., 2006). Ver `classificacao_koppen_grade`.

    Parâmetros
    ----------
    T : lista ou array de 12 valores
        Temperatura média mensal do ar (jan-dez), em °C.
    P : lista ou array de 12 valores
        Precipitação mensal (jan-dez), em mm/mês.
    lat : float
        Latitude, em graus (negativa no hemisfério sul).

    Retorna
    -------
    dict
        'id' (int) e 'classe' (str, ex.: 'Cfa').
    """
    _um_local(T, P)
    return _escalar(classificacao_koppen_grade(T, P, lat))


# ---------------------------------------------------------------------------
# Thornthwaite (1948)
# ---------------------------------------------------------------------------
TH_UMIDADE = {1: "A", 2: "B4", 3: "B3", 4: "B2", 5: "B1", 6: "C2", 7: "C1", 8: "D", 9: "E"}
TH_UMIDADE_NOMES = {1: "Superúmido", 2: "Úmido", 3: "Úmido", 4: "Úmido", 5: "Úmido",
                    6: "Subúmido", 7: "Subúmido seco", 8: "Semiárido", 9: "Árido"}
TH_SUBTIPO = {1: "r", 2: "s", 3: "w", 4: "s2", 5: "w2",   # úmidos (Im >= 0), pela deficiência
              6: "d", 7: "s", 8: "w", 9: "s2", 10: "w2"}  # secos (Im < 0), pelo excedente
TH_TERMICA = {1: "A'", 2: "B'4", 3: "B'3", 4: "B'2", 5: "B'1", 6: "C'2", 7: "C'1", 8: "D'", 9: "E'"}
TH_TERMICA_NOMES = {1: "Megatérmico", 2: "Mesotérmico", 3: "Mesotérmico", 4: "Mesotérmico",
                    5: "Mesotérmico", 6: "Microtérmico", 7: "Microtérmico", 8: "Tundra",
                    9: "Gelo perpétuo"}
TH_CONCENTRACAO = {1: "a'", 2: "b'4", 3: "b'3", 4: "b'2", 5: "b'1", 6: "c'2", 7: "c'1", 8: "d'"}


def classificacao_thornthwaite_grade(P, ETP, lat, CAD=100.0):
    """
    Classificação climática de Thornthwaite (1948) para muitos locais, a
    partir do balanço hídrico climatológico de Thornthwaite & Mather (1955).

    Tabelas de Aparecido et al. (2016), com três correções pelo original de
    Thornthwaite (1948): limite B'3/B'2 em 855 mm (o artigo traz 885); nos
    climas secos, s2/w2 a partir de Ih >= 20 (o artigo traz 33,3) e excedente
    maior no verão = w, no inverno = s (o artigo traz invertido).

    Índices (%): Ih = 100 EXC/ETP, Ia = 100 DEF/ETP, Im = Ih - 0,6 Ia.
    Concentração estival: ETP do verão (1/3 dez + jan + fev + 2/3 mar, no
    hemisfério sul) em % da ETP anual.

    Parâmetros
    ----------
    P, ETP : array_like, formato (12, ...)
        Precipitação e evapotranspiração potencial mensais, em mm/mês.
    lat : float ou array_like
        Latitude, em graus (define o hemisfério para verão e inverno).
    CAD : float ou array_like, opcional
        Capacidade de água disponível no solo, em mm. Padrão 100.

    Retorna
    -------
    dict de numpy.ndarray
        'classe' (ex.: "B1rA'a'"), 'umidade', 'subtipo', 'termica',
        'concentracao' (códigos inteiros, 0 = sem dado), 'Ih', 'Ia', 'Im',
        'ETP_anual', 'DEF_anual', 'EXC_anual' (mm) e 'ETP_verao_pct'.
    """
    P = _mensal(P, "P")
    ETP = _mensal(ETP, "ETP")
    forma = P.shape[1:]
    norte = _espacial(lat, forma) > 0
    CAD = _espacial(CAD, forma)
    valido = (np.all(np.isfinite(P), axis=0) & np.all(np.isfinite(ETP), axis=0)
              & np.isfinite(CAD) & (CAD > 0))

    bh = balanco_hidrico_climatologico_grade(P, ETP, np.where(valido, CAD, 1.0))
    etp_anual = ETP.sum(0)
    etp_seguro = np.maximum(etp_anual, 1e-9)
    def_anual, exc_anual = bh['DEF'].sum(0), bh['EXC'].sum(0)
    ih = 100 * exc_anual / etp_seguro
    ia = 100 * def_anual / etp_seguro
    im = ih - 0.6 * ia
    petr = 100 * _soma_estacao(ETP, 'verao', norte) / etp_seguro

    umidade = _faixas(im, [-40, -20, 0, 20, 40, 60, 80, 100], [9, 8, 7, 6, 5, 4, 3, 2, 1])

    def_v, def_i = _soma_estacao(bh['DEF'], 'verao', norte), _soma_estacao(bh['DEF'], 'inverno', norte)
    exc_v, exc_i = _soma_estacao(bh['EXC'], 'verao', norte), _soma_estacao(bh['EXC'], 'inverno', norte)
    umido, seco = im >= 0, im < 0
    def_verao, exc_verao = def_v > def_i, exc_v > exc_i
    subtipo = np.zeros(forma, dtype=np.uint8)
    for cond, c in [
        (umido & (ia < 16.7), 1),
        (umido & (ia >= 16.7) & (ia < 33.3) & def_verao, 2),
        (umido & (ia >= 16.7) & (ia < 33.3) & ~def_verao, 3),
        (umido & (ia >= 33.3) & def_verao, 4),
        (umido & (ia >= 33.3) & ~def_verao, 5),
        (seco & (ih < 10), 6),
        (seco & (ih >= 10) & (ih < 20) & ~exc_verao, 7),
        (seco & (ih >= 10) & (ih < 20) & exc_verao, 8),
        (seco & (ih >= 20) & ~exc_verao, 9),
        (seco & (ih >= 20) & exc_verao, 10),
    ]:
        subtipo[cond] = c

    termica = _faixas(etp_anual, [142, 285, 427, 570, 712, 855, 997, 1140], [9, 8, 7, 6, 5, 4, 3, 2, 1])
    concentracao = _faixas(petr, [48, 51.9, 56.3, 61.6, 68, 76.3, 88], [1, 2, 3, 4, 5, 6, 7, 8])

    for c in (umidade, subtipo, termica, concentracao):
        c[~valido] = 0
    classe = (_rotulos(umidade, TH_UMIDADE) + _rotulos(subtipo, TH_SUBTIPO)
              + _rotulos(termica, TH_TERMICA) + _rotulos(concentracao, TH_CONCENTRACAO))
    classe = np.where(valido, classe, '').astype(str)

    nan = np.where(valido, 1.0, np.nan)
    return {
        'classe': classe,
        'umidade': umidade, 'subtipo': subtipo, 'termica': termica, 'concentracao': concentracao,
        'Ih': ih * nan, 'Ia': ia * nan, 'Im': im * nan,
        'ETP_anual': etp_anual * nan, 'DEF_anual': def_anual * nan, 'EXC_anual': exc_anual * nan,
        'ETP_verao_pct': petr * nan,
    }


def classificacao_thornthwaite(P, ETP, lat, CAD=100.0):
    """
    Classificação climática de Thornthwaite (1948) de um local. Ver
    `classificacao_thornthwaite_grade`.

    Parâmetros
    ----------
    P, ETP : lista ou array de 12 valores
        Precipitação e evapotranspiração potencial mensais (jan-dez), em mm/mês.
    lat : float
        Latitude, em graus (negativa no hemisfério sul).
    CAD : float, opcional
        Capacidade de água disponível no solo, em mm. Padrão 100.

    Retorna
    -------
    dict
        'classe' (ex.: "B1rA'a'"), 'descricao', códigos e índices (Ih, Ia, Im...).
    """
    _um_local(P, ETP)
    saida = _escalar(classificacao_thornthwaite_grade(P, ETP, lat, CAD))
    if saida['classe']:
        saida['descricao'] = (f"{TH_UMIDADE_NOMES[saida['umidade']]}, "
                              f"{TH_TERMICA_NOMES[saida['termica']].lower()}")
    else:
        saida['descricao'] = ''
    return saida


# ---------------------------------------------------------------------------
# Camargo (1991) modificado por Maluf (2000)
# ---------------------------------------------------------------------------
CAMARGO_TERMICA = {1: 'GL', 2: 'FR', 3: 'CO', 4: 'TE', 5: 'STE', 6: 'ST', 7: 'TR', 8: 'EQ'}
CAMARGO_TERMICA_NOMES = {1: 'Glacial', 2: 'Frio', 3: 'Frio moderado', 4: 'Temperado',
                         5: 'Subtemperado', 6: 'Subtropical', 7: 'Tropical', 8: 'Equatorial'}
CAMARGO_HIDRICA = {1: 'DE', 2: 'AR', 3: 'SE', 4: 'MO', 5: 'SB', 6: 'UM', 7: 'PU', 8: 'SU'}
CAMARGO_HIDRICA_NOMES = {1: 'Desértico', 2: 'Árido', 3: 'Seco', 4: 'Monçônico', 5: 'Subúmido',
                         6: 'Úmido', 7: 'Superúmido', 8: 'Extremamente úmido'}
CAMARGO_ESTACAO_SECA = {1: 'v', 2: 'o', 3: 'i', 4: 'p'}  # verão, outono, inverno, primavera

_ZERO_MM = 0.05  # DEF e EXC abaixo disso contam como zero (arredondam para 0,0 mm)


def classificacao_camargo_grade(T, P, ETP, lat, CAD=100.0, tabela_termica='coerente'):
    """
    Classificação climática de Camargo (1991) modificada por Maluf (2000)
    para muitos locais (Tabelas 6-8 de Aparecido et al., 2016).

    Classe térmica pela temperatura média anual (Ty) e do mês mais frio
    (Tcold); classe hídrica pela deficiência (DEF) e pelo excedente (EXC)
    anuais do balanço hídrico climatológico; letra da estação com maior
    deficiência (v, o, i, p) nas classes SE, MO, SB e UM.

    A Tabela 6 do artigo é ambígua para Tcold > 20 °C. Com
    `tabela_termica='coerente'` (padrão), Tcold só subdivide a faixa
    18 < Ty <= 22 (STE <= 13 < ST <= 20 < TR) e Ty > 25 é sempre EQ. Com
    'literal', Tcold > 20 °C leva a TR em qualquer caso (como no script do
    GEE do repositório climas_brasil), e a classe EQ praticamente não ocorre.

    A Tabela 7 não cobre 0 < DEF <= 150 mm com EXC = 0 nem DEF = 0 com
    EXC <= 200 mm; esses casos entram em SB (0 <= DEF <= 150 e 0 <= EXC <= 200).

    Parâmetros
    ----------
    T : array_like, formato (12, ...)
        Temperatura média mensal do ar, em °C.
    P, ETP : array_like, formato (12, ...)
        Precipitação e evapotranspiração potencial mensais, em mm/mês.
    lat : float ou array_like
        Latitude, em graus (define o hemisfério das estações).
    CAD : float ou array_like, opcional
        Capacidade de água disponível no solo, em mm. Padrão 100.
    tabela_termica : {'coerente', 'literal'}, opcional
        Leitura da Tabela 6 (ver acima). Padrão 'coerente'.

    Retorna
    -------
    dict de numpy.ndarray
        'classe' (ex.: 'ST-UMi'), 'termica', 'hidrica', 'estacao_seca'
        (códigos inteiros, 0 = sem dado/sem letra), 'T_anual', 'T_mes_frio',
        'DEF_anual' e 'EXC_anual'.
    """
    if tabela_termica not in ('coerente', 'literal'):
        raise ValueError("tabela_termica deve ser 'coerente' ou 'literal'.")
    T = _mensal(T, "T")
    P = _mensal(P, "P")
    ETP = _mensal(ETP, "ETP")
    forma = T.shape[1:]
    norte = _espacial(lat, forma) > 0
    CAD = _espacial(CAD, forma)
    valido = (np.all(np.isfinite(T), axis=0) & np.all(np.isfinite(P), axis=0)
              & np.all(np.isfinite(ETP), axis=0) & np.isfinite(CAD) & (CAD > 0))

    ty = T.mean(0)
    tcold = T.min(0)
    bh = balanco_hidrico_climatologico_grade(P, ETP, np.where(valido, CAD, 1.0))
    def_anual, exc_anual = bh['DEF'].sum(0), bh['EXC'].sum(0)
    DEF = np.where(def_anual < _ZERO_MM, 0.0, def_anual)
    EXC = np.where(exc_anual < _ZERO_MM, 0.0, exc_anual)

    termica = np.zeros(forma, dtype=np.uint8)
    for cond, c in [
        (ty <= 3, 1),
        ((ty > 3) & (ty <= 7), 2),
        ((ty > 7) & (ty <= 12), 3),
        ((ty > 12) & (ty <= 18), 4),
        ((ty > 18) & (ty <= 22) & (tcold <= 13), 5),
        ((ty > 18) & (ty <= 22) & (tcold > 13) & (tcold <= 20), 6),
        ((ty > 18) & (ty <= 22) & (tcold > 20), 7),
        ((ty > 22) & (ty <= 25), 7),
        (ty > 25, 8),
    ]:
        termica[cond] = c
    if tabela_termica == 'literal':
        termica[tcold > 20] = 7

    hidrica = np.zeros(forma, dtype=np.uint8)
    for cond, c in [
        ((DEF > 800) & (EXC == 0), 1),
        ((DEF > 150) & (DEF <= 800) & (EXC == 0), 2),
        ((DEF > 150) & (EXC > 0) & (EXC <= 200), 3),
        ((DEF > 150) & (EXC > 200), 4),
        ((DEF <= 150) & (EXC <= 200), 5),
        ((DEF > 0) & (DEF <= 150) & (EXC > 200), 6),
        ((DEF == 0) & (EXC > 200) & (EXC <= 1000), 7),
        ((DEF == 0) & (EXC > 1000), 8),
    ]:
        hidrica[cond] = c

    def_estacoes = np.stack([_soma_estacao(bh['DEF'], e, norte)
                             for e in ('verao', 'outono', 'inverno', 'primavera')])
    com_letra = np.isin(hidrica, [3, 4, 5, 6]) & (def_estacoes.max(0) >= _ZERO_MM)
    estacao_seca = np.where(com_letra, def_estacoes.argmax(0) + 1, 0).astype(np.uint8)

    for c in (termica, hidrica, estacao_seca):
        c[~valido] = 0
    classe = (_rotulos(termica, CAMARGO_TERMICA) + '-' + _rotulos(hidrica, CAMARGO_HIDRICA)
              + _rotulos(estacao_seca, CAMARGO_ESTACAO_SECA))
    classe = np.where(valido, classe, '').astype(str)

    nan = np.where(valido, 1.0, np.nan)
    return {
        'classe': classe, 'termica': termica, 'hidrica': hidrica, 'estacao_seca': estacao_seca,
        'T_anual': ty * nan, 'T_mes_frio': tcold * nan,
        'DEF_anual': def_anual * nan, 'EXC_anual': exc_anual * nan,
    }


def classificacao_camargo(T, P, ETP, lat, CAD=100.0, tabela_termica='coerente'):
    """
    Classificação climática de Camargo (1991) mod. Maluf (2000) de um local.
    Ver `classificacao_camargo_grade`.

    Parâmetros
    ----------
    T : lista ou array de 12 valores
        Temperatura média mensal do ar (jan-dez), em °C.
    P, ETP : lista ou array de 12 valores
        Precipitação e evapotranspiração potencial mensais, em mm/mês.
    lat : float
        Latitude, em graus (negativa no hemisfério sul).
    CAD : float, opcional
        Capacidade de água disponível no solo, em mm. Padrão 100.
    tabela_termica : {'coerente', 'literal'}, opcional
        Leitura da Tabela 6 de Aparecido et al. (2016). Padrão 'coerente'.

    Retorna
    -------
    dict
        'classe' (ex.: 'ST-UMi'), 'descricao', códigos e variáveis anuais.
    """
    _um_local(T, P, ETP)
    saida = _escalar(classificacao_camargo_grade(T, P, ETP, lat, CAD, tabela_termica))
    if saida['classe']:
        saida['descricao'] = (f"{CAMARGO_TERMICA_NOMES[saida['termica']]} "
                              f"{CAMARGO_HIDRICA_NOMES[saida['hidrica']].lower()}")
    else:
        saida['descricao'] = ''
    return saida


# ---------------------------------------------------------------------------
# Zonas de vida de Holdridge (38 zonas; Leemans, 1990; Jungkunst et al., 2021)
# ---------------------------------------------------------------------------
HOLDRIDGE_ZONAS = {
    1: "Gelo polar", 2: "Deserto polar",
    3: "Tundra seca subpolar", 4: "Tundra úmida subpolar", 5: "Tundra muito úmida subpolar",
    6: "Tundra pluvial subpolar",
    7: "Deserto boreal", 8: "Arbustal seco boreal", 9: "Floresta úmida boreal",
    10: "Floresta muito úmida boreal", 11: "Floresta pluvial boreal",
    12: "Deserto temperado frio", 13: "Arbustal desértico temperado frio",
    14: "Estepe temperada fria", 15: "Floresta úmida temperada fria",
    16: "Floresta muito úmida temperada fria", 17: "Floresta pluvial temperada fria",
    18: "Deserto temperado quente", 19: "Arbustal desértico temperado quente",
    20: "Estepe espinhosa temperada quente", 21: "Floresta seca temperada quente",
    22: "Floresta úmida temperada quente", 23: "Floresta muito úmida temperada quente",
    24: "Floresta pluvial temperada quente",
    25: "Deserto subtropical", 26: "Arbustal desértico subtropical",
    27: "Estepe espinhosa subtropical", 28: "Floresta seca subtropical",
    29: "Floresta úmida subtropical", 30: "Floresta muito úmida subtropical",
    31: "Floresta pluvial subtropical",
    32: "Deserto tropical", 33: "Arbustal desértico tropical", 34: "Estepe espinhosa tropical",
    35: "Floresta muito seca tropical", 36: "Floresta seca tropical", 37: "Floresta úmida tropical",
    38: "Floresta muito úmida tropical",
}

# Nomes originais (Jungkunst et al., 2021, Tab. 1).
HOLDRIDGE_ZONAS_EN = {
    1: "Polar ice", 2: "Polar desert",
    3: "Subpolar dry tundra", 4: "Subpolar moist tundra", 5: "Subpolar wet tundra",
    6: "Subpolar rain tundra",
    7: "Boreal desert", 8: "Boreal dry bush", 9: "Boreal moist forest", 10: "Boreal wet forest",
    11: "Boreal rain forest",
    12: "Cool temperate desert", 13: "Cool temperate desert bush", 14: "Cool temperate steppe",
    15: "Cool temperate moist forest", 16: "Cool temperate wet forest",
    17: "Cool temperate rain forest",
    18: "Warm temperate desert", 19: "Warm temperate desert bush", 20: "Warm temperate thorn steppe",
    21: "Warm temperate dry forest", 22: "Warm temperate moist forest",
    23: "Warm temperate wet forest", 24: "Warm temperate rain forest",
    25: "Subtropical desert", 26: "Subtropical desert bush", 27: "Subtropical thorn steppe",
    28: "Subtropical dry forest", 29: "Subtropical moist forest", 30: "Subtropical wet forest",
    31: "Subtropical rain forest",
    32: "Tropical desert", 33: "Tropical desert bush", 34: "Tropical thorn steppe",
    35: "Tropical very dry forest", 36: "Tropical dry forest", 37: "Tropical moist forest",
    38: "Tropical wet forest",
}

# (classe de biotemperatura, classe de ETP/P mínima, máxima, zona). Cada faixa térmica começa
# numa razão ETP/P diferente (Jungkunst et al., 2021, Fig. 1); valores além das pontas entram
# na zona extrema da faixa.
_HOLDRIDGE_TABELA = [
    (2, 1, 6, 3), (2, 7, 7, 4), (2, 8, 8, 5), (2, 9, 10, 6),
    (3, 1, 5, 7), (3, 6, 6, 8), (3, 7, 7, 9), (3, 8, 8, 10), (3, 9, 10, 11),
    (4, 1, 4, 12), (4, 5, 5, 13), (4, 6, 6, 14), (4, 7, 7, 15), (4, 8, 8, 16), (4, 9, 10, 17),
    (5, 1, 3, 18), (5, 4, 4, 19), (5, 5, 5, 20), (5, 6, 6, 21), (5, 7, 7, 22), (5, 8, 8, 23),
    (5, 9, 10, 24),
    (6, 1, 3, 25), (6, 4, 4, 26), (6, 5, 5, 27), (6, 6, 6, 28), (6, 7, 7, 29), (6, 8, 8, 30),
    (6, 9, 10, 31),
    (7, 1, 2, 32), (7, 3, 3, 33), (7, 4, 4, 34), (7, 5, 5, 35), (7, 6, 6, 36), (7, 7, 7, 37),
    (7, 8, 10, 38),
]

FATOR_ETP_HOLDRIDGE = 58.93  # mm/ano por °C de biotemperatura (Holdridge, 1967)


def classificacao_holdridge_grade(T, P, lat, ETP=None, limiar_correcao=24.0):
    """
    Zonas de vida de Holdridge (38 zonas) para muitos locais, pela
    biotemperatura e pela razão de evapotranspiração potencial (ETP/P).

    Biotemperatura: média das temperaturas mensais corrigidas pela
    latitude, t - 3|lat|/100 (t - 24)^2 nos meses com t > `limiar_correcao`,
    e limitadas a [0, 30] °C.

    Parâmetros
    ----------
    T : array_like, formato (12, ...)
        Temperatura média mensal do ar, em °C.
    P : array_like, formato (12, ...)
        Precipitação mensal, em mm/mês.
    lat : float ou array_like
        Latitude, em graus.
    ETP : array_like, formato (12, ...), opcional
        ETP mensal, em mm/mês (ex.: Penman-Monteith). Se None (padrão), usa
        a ETP de Holdridge: 58,93 x biotemperatura (mm/ano).
    limiar_correcao : float ou None, opcional
        Temperatura mensal (°C) acima da qual se aplica a correção por
        latitude. Padrão 24 (Holdridge, 1967). None aplica a todos os meses.

    Retorna
    -------
    dict de numpy.ndarray
        'id' (1-38, ver HOLDRIDGE_ZONAS; 0 = sem dado), 'classe' (nome da
        zona em português), 'classe_en' (nome original em inglês),
        'biotemperatura' (°C), 'P_anual' e 'ETP_anual' (mm) e
        'razao_ETP' (ETP/P).
    """
    T = _mensal(T, "T")
    P = _mensal(P, "P")
    forma = T.shape[1:]
    alat = np.abs(_espacial(lat, forma))
    valido = np.all(np.isfinite(T), axis=0) & np.all(np.isfinite(P), axis=0)

    correcao = 3 * alat / 100 * (T - 24) ** 2
    if limiar_correcao is not None:
        correcao = np.where(T > limiar_correcao, correcao, 0.0)
    biotemp = np.clip(T - correcao, 0, 30).mean(0)

    p_anual = P.sum(0)
    if ETP is None:
        etp_anual = FATOR_ETP_HOLDRIDGE * biotemp
    else:
        ETP = _mensal(ETP, "ETP")
        etp_anual = ETP.sum(0)
        valido &= np.all(np.isfinite(ETP), axis=0)
    razao = etp_anual / np.maximum(p_anual, 1.0)

    t_cls = _faixas(biotemp, [1.5, 3, 6, 12, 18, 24], [1, 2, 3, 4, 5, 6, 7])
    u_cls = _faixas(razao, [0.125, 0.25, 0.5, 1, 2, 4, 8, 16, 32], [10, 9, 8, 7, 6, 5, 4, 3, 2, 1])

    zona = np.zeros(forma, dtype=np.uint8)
    zona[(t_cls == 1) & (biotemp < 1)] = 1
    zona[(t_cls == 1) & (biotemp >= 1)] = 2
    for tc, umin, umax, z in _HOLDRIDGE_TABELA:
        zona[(t_cls == tc) & (u_cls >= umin) & (u_cls <= umax)] = z
    zona[~valido] = 0

    nan = np.where(valido, 1.0, np.nan)
    return {
        'id': zona, 'classe': _rotulos(zona, HOLDRIDGE_ZONAS).astype(str),
        'classe_en': _rotulos(zona, HOLDRIDGE_ZONAS_EN).astype(str),
        'biotemperatura': biotemp * nan, 'P_anual': p_anual * nan,
        'ETP_anual': etp_anual * nan, 'razao_ETP': razao * nan,
    }


def classificacao_holdridge(T, P, lat, ETP=None, limiar_correcao=24.0):
    """
    Zona de vida de Holdridge de um local. Ver `classificacao_holdridge_grade`.

    Parâmetros
    ----------
    T : lista ou array de 12 valores
        Temperatura média mensal do ar (jan-dez), em °C.
    P : lista ou array de 12 valores
        Precipitação mensal, em mm/mês.
    lat : float
        Latitude, em graus.
    ETP : lista ou array de 12 valores, opcional
        ETP mensal, em mm/mês. Se None, usa 58,93 x biotemperatura.
    limiar_correcao : float ou None, opcional
        Limiar (°C) da correção por latitude. Padrão 24.

    Retorna
    -------
    dict
        'id', 'classe' (nome da zona), 'classe_en', 'biotemperatura', 'P_anual',
        'ETP_anual' e 'razao_ETP'.
    """
    _um_local(T, P, ETP)
    return _escalar(classificacao_holdridge_grade(T, P, lat, ETP, limiar_correcao))
