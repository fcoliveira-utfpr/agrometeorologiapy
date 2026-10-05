"""Capítulo 12 — Produtividade potencial e atingível (Zona Agroecológica da FAO).

Modelo de Doorenbos & Kassam (1979), com as correções de temperatura de
Barbieri & Tuon (1992) (em Pereira et al., 2002) e a razão de insolação
obtida pela inversão da equação de Angström-Prescott.
"""

from calendar import monthrange

import numpy as np
import pandas as pd

from ._trig import cosd
from .radiacao import angulo_horario_nascer, declinacao_solar, fator_correcao_distancia, irradiancia_extraterrestre

__all__ = [
    "razao_insolacao",
    "cTn",
    "cTc",
    "PPBp",
    "CIAF",
    "CR",
    "produtividade_potencial",
    "produtividade_atingivel",
]


def _como_entrada(valor, ref):
    """Devolve `valor` no mesmo tipo de `ref`: float para escalar, Series
    com o mesmo índice para Series, numpy.ndarray nos demais casos."""
    if isinstance(ref, pd.Series):
        return pd.Series(np.asarray(valor, dtype=float), index=ref.index)
    if np.ndim(ref) == 0:
        return float(valor)
    return np.asarray(valor, dtype=float)


def _rota(rota):
    r = str(rota).upper()
    if r not in ("C3", "C4"):
        raise ValueError("rota deve ser 'C3' ou 'C4'.")
    return r


def razao_insolacao(Qg, Qo, lat, a=None, b=0.52):
    """
    Razão de insolação (n/N) pela inversão da equação de Angström-Prescott,
    para quando não há heliógrafo (ex.: bases em grade como o BR-DWGD).

    n/N = (Qg/Qo - a) / b, limitada entre 0 e 1.

    Parâmetros
    ----------
    Qg : float ou array_like
        Irradiância solar global, em MJ/m² dia.
    Qo : float ou array_like
        Irradiância solar extraterrestre, em MJ/m² dia.
    lat : float
        Latitude do local, em graus (negativa no hemisfério sul).
    a : float, opcional
        Coeficiente linear de Angström-Prescott. Padrão None, que usa
        a = 0,29 cos(lat) (Glover-McCulloch), como em `Qg_angstrom`.
    b : float, opcional
        Coeficiente angular de Angström-Prescott. Padrão 0,52.

    Retorna
    -------
    nN : float ou array_like
        Razão de insolação (adimensional, 0 a 1).

    Exemplos
    --------
    >>> round(razao_insolacao(22.5, 42.96, -24.86), 3)
    0.501
    """
    if a is None:
        a = 0.29 * cosd(lat)
    nN = np.clip((np.asarray(Qg, dtype=float) / np.asarray(Qo, dtype=float) - a) / b, 0.0, 1.0)
    return _como_entrada(nN, Qg)


def cTn(T, rota, T_limiar=16.5):
    """
    Correção de temperatura da fotossíntese para o céu nublado (cTn),
    pelos polinômios de Barbieri & Tuon (1992), limitada a zero.

    C3: cTn = -0,0425 + 0,035 T + 0,00325 T² - 0,0000925 T³
    C4: cTn = -1,064 + 0,173 T - 0,0029 T²       (T >= 16,5 °C)
        cTn = -4,16 + 0,4325 T - 0,00725 T²      (T < 16,5 °C)

    Parâmetros
    ----------
    T : float ou array_like
        Temperatura média do ar no período, em °C.
    rota : str
        Rota fotossintética da cultura: 'C3' (soja, trigo, feijão,
        girassol, arroz...) ou 'C4' (milho, sorgo, cana-de-açúcar...).
    T_limiar : float, opcional
        Temperatura que separa os dois polinômios das culturas C4, em °C.
        Padrão 16,5. Ignorado para C3.

    Retorna
    -------
    cTn : float ou array_like
        Correção de temperatura para céu nublado (adimensional, >= 0).

    Exemplos
    --------
    >>> round(cTn(25, 'C3'), 3), round(cTn(25, 'C4'), 3)
    (1.418, 1.448)
    """
    rota = _rota(rota)
    t = np.asarray(T, dtype=float)
    if rota == "C4":
        c = np.where(t >= T_limiar,
                     -1.064 + 0.173 * t - 0.0029 * t ** 2,
                     -4.16 + 0.4325 * t - 0.00725 * t ** 2)
    else:
        c = -0.0425 + 0.035 * t + 0.00325 * t ** 2 - 0.0000925 * t ** 3
    return _como_entrada(np.clip(c, 0.0, None), T)


def cTc(T, rota, T_limiar=16.5):
    """
    Correção de temperatura da fotossíntese para o céu claro (cTc),
    pelos polinômios de Barbieri & Tuon (1992), limitada a zero.

    C3: cTc = 0,583 + 0,014 T + 0,0013 T² - 0,000037 T³
    C4: cTc = -4,16 + 0,4325 T - 0,00725 T²      (T >= 16,5 °C)
        cTc = -9,32 + 0,865 T - 0,0145 T²        (T < 16,5 °C)

    Parâmetros
    ----------
    T : float ou array_like
        Temperatura média do ar no período, em °C.
    rota : str
        Rota fotossintética da cultura: 'C3' ou 'C4'.
    T_limiar : float, opcional
        Temperatura que separa os dois polinômios das culturas C4, em °C.
        Padrão 16,5. Ignorado para C3.

    Retorna
    -------
    cTc : float ou array_like
        Correção de temperatura para céu claro (adimensional, >= 0).

    Exemplos
    --------
    >>> round(cTc(25, 'C3'), 3), round(cTc(25, 'C4'), 3)
    (1.167, 2.121)
    """
    rota = _rota(rota)
    t = np.asarray(T, dtype=float)
    if rota == "C4":
        c = np.where(t >= T_limiar,
                     -4.16 + 0.4325 * t - 0.00725 * t ** 2,
                     -9.32 + 0.865 * t - 0.0145 * t ** 2)
    else:
        c = 0.583 + 0.014 * t + 0.0013 * t ** 2 - 0.000037 * t ** 3
    return _como_entrada(np.clip(c, 0.0, None), T)


def PPBp(Qo, nN, cTc, cTn, a_c=107.2, b_c=8.604, a_n=31.7, b_n=5.234):
    """
    Produtividade potencial bruta padrão (PPBp) da cultura padrão
    hipotética (IAF = 5, sem restrição de água e nutrientes), soma das
    frações de céu claro e de céu nublado do dia.

    PPBc = (a_c + b_c Qo) cTc (n/N)
    PPBn = (a_n + b_n Qo) cTn (1 - n/N)
    PPBp = PPBc + PPBn

    Os coeficientes padrão são os de Doorenbos & Kassam convertidos para
    Qo em MJ/m² dia (0,36 e 0,219 com Qo em cal/cm² dia).

    Parâmetros
    ----------
    Qo : float ou array_like
        Irradiância solar extraterrestre, em MJ/m² dia.
    nN : float ou array_like
        Razão de insolação n/N (0 a 1). Ver `razao_insolacao`.
    cTc, cTn : float ou array_like
        Correções de temperatura para céu claro e nublado. Ver `cTc` e `cTn`.
    a_c, b_c : float, opcional
        Coeficientes da fração de céu claro. Padrão 107,2 e 8,604.
    a_n, b_n : float, opcional
        Coeficientes da fração de céu nublado. Padrão 31,7 e 5,234.

    Retorna
    -------
    PPBp : float ou array_like
        Produtividade potencial bruta padrão, em kg MS/ha dia.

    Exemplos
    --------
    >>> round(PPBp(42.96, 0.501, cTc=2.210, cTn=1.484), 0)
    718.0
    """
    PPBc = (a_c + b_c * np.asarray(Qo, dtype=float)) * cTc * nN
    PPBn = (a_n + b_n * np.asarray(Qo, dtype=float)) * cTn * (1 - np.asarray(nN, dtype=float))
    return _como_entrada(PPBc + PPBn, Qo)


def CIAF(IAF, IAF_padrao=5.0):
    """
    Correção do índice de área foliar (CIAF): ajusta a cultura padrão
    (IAF = 5) ao IAF máximo da cultura.

    CIAF = 0,0093 + 0,185 IAF - 0,0175 IAF²   (IAF < 5)
    CIAF = 0,5                                (IAF >= 5)

    Parâmetros
    ----------
    IAF : float
        Índice de área foliar máximo da cultura (m²/m²).
    IAF_padrao : float, opcional
        IAF a partir do qual CIAF = 0,5. Padrão 5,0.

    Retorna
    -------
    CIAF : float
        Correção do índice de área foliar (adimensional).

    Exemplos
    --------
    >>> round(CIAF(4.5), 3), CIAF(5)
    (0.487, 0.5)
    """
    i = np.asarray(IAF, dtype=float)
    c = np.where(i >= IAF_padrao, 0.5, 0.0093 + 0.185 * i - 0.0175 * i ** 2)
    return _como_entrada(c, IAF)


def CR(T, T_limiar=20.0, CR_frio=0.6, CR_quente=0.5):
    """
    Correção para a respiração de manutenção (CR), conforme a temperatura
    média do período.

    CR = 0,6 se T < 20 °C;  CR = 0,5 se T >= 20 °C.

    Parâmetros
    ----------
    T : float ou array_like
        Temperatura média do ar no período, em °C.
    T_limiar : float, opcional
        Temperatura de troca do CR, em °C. Padrão 20,0.
    CR_frio, CR_quente : float, opcional
        CR abaixo e a partir de `T_limiar`. Padrão 0,6 e 0,5.

    Retorna
    -------
    CR : float ou array_like
        Correção para a respiração (adimensional).

    Exemplos
    --------
    >>> CR(18.0), CR(26.5)
    (0.6, 0.5)
    """
    c = np.where(np.asarray(T, dtype=float) < T_limiar, CR_frio, CR_quente)
    return _como_entrada(c, T)


def _inicio_periodos(df, intervalo, ano):
    """Data de início de cada período (pandas.DatetimeIndex) ou None, se o
    df não traz datas. Aceita a coluna 'data', um DatetimeIndex ou as
    colunas 'dia' e 'mes' (+ 'ano' opcional; sem ela, o ano começa em `ano`
    e avança quando o mês diminui)."""
    if 'data' in df.columns:
        datas = pd.DatetimeIndex(pd.to_datetime(df['data']))
    elif isinstance(df.index, pd.DatetimeIndex):
        datas = df.index
    elif {'dia', 'mes'} <= set(df.columns):
        meses = df['mes'].astype(int).to_numpy()
        if 'ano' in df.columns:
            anos = df['ano'].astype(int).to_numpy()
        else:
            anos = ano + np.concatenate([[0], np.cumsum(np.diff(meses) < 0)])
        datas = pd.DatetimeIndex([pd.Timestamp(int(a), int(m), int(d))
                                  for a, m, d in zip(anos, meses, df['dia'].astype(int))])
    else:
        return None

    # Normaliza para o início do período do calendário.
    if intervalo == 'dec':
        dia_ini = np.minimum((datas.day - 1) // 10, 2) * 10 + 1
        datas = pd.DatetimeIndex([pd.Timestamp(d.year, d.month, int(di)) for d, di in zip(datas, dia_ini)])
    elif intervalo == 'M':
        datas = pd.DatetimeIndex([pd.Timestamp(d.year, d.month, 1) for d in datas])
    return datas.normalize()


def _dias_periodo(inicio, intervalo):
    """Número de dias do período do calendário que começa em `inicio`."""
    if intervalo == 'd':
        return 1
    dias_mes = monthrange(inicio.year, inicio.month)[1]
    if intervalo == 'M':
        return dias_mes
    return 10 if inicio.day <= 11 else dias_mes - 20


def produtividade_potencial(df, lat, rota, IAF=5.0, Cc=0.35, U=13.0, intervalo='dec',
                            semeadura=None, ciclo=None, ano=2023, a=None, b=0.52,
                            T_limiar_CR=20.0, CR_frio=0.6, CR_quente=0.5):
    """
    Produtividade potencial (PPf) pelo Modelo da Zona Agroecológica da FAO
    (Doorenbos & Kassam, 1979), acumulada período a período ao longo do
    ciclo da cultura.

    Em cada período (dia, decêndio ou mês):

    1. n/N pela inversão de Angström-Prescott (se não for informada);
    2. cTn e cTc da rota C3 ou C4 com a temperatura média do período;
    3. PPBp = PPBc + PPBn (kg MS/ha dia);
    4. CR conforme a temperatura do período;
    5. PPf_periodo = PPBp x ND x CIAF x CR x Cc / (1 - 0,01 U),
       com ND = dias do período em que a cultura estava no campo.

    PPf = soma de PPf_periodo no ciclo (coluna 'PPf_acum').

    Parâmetros
    ----------
    df : pandas.DataFrame
        Uma linha por período, em ordem cronológica. Colunas:

        - 'Tmed' (°C), obrigatória: temperatura média do período.
        - Radiação: 'nN' (razão de insolação), ou 'n' e 'N' (horas), ou
          'Qg' (MJ/m² dia, média diária do período), nessa prioridade.
        - 'Qo' (MJ/m² dia, média diária do período), opcional: se faltar,
          é calculada a partir de `lat` e das datas.
        - Datas: coluna 'data', DatetimeIndex ou colunas 'dia' e 'mes'
          (+ 'ano' opcional), com o início de cada período. Necessárias
          para calcular Qo, ND e para usar `semeadura`.
        - 'ND', opcional: número de dias do período a considerar. Se
          faltar, vem do calendário (decêndio: 10, 10 e 8 a 11 dias).
    lat : float
        Latitude do local, em graus (negativa no hemisfério sul).
    rota : str
        Rota fotossintética da cultura: 'C3' ou 'C4'.
    IAF : float, opcional
        Índice de área foliar máximo da cultura. Padrão 5,0.
    Cc : float, opcional
        Coeficiente (índice) de colheita. Padrão 0,35.
    U : float, opcional
        Umidade do produto colhido, em %. Padrão 13,0.
    intervalo : str, opcional
        Escala de tempo das linhas do df: 'd' (diário), 'dec' (decendial
        do calendário: 1-10, 11-20, 21-fim do mês) ou 'M' (mensal).
        Padrão 'dec'.
    semeadura : str ou datetime, opcional
        Data de semeadura ('AAAA-MM-DD'). Com `ciclo`, só os dias de
        semeadura a semeadura + ciclo - 1 entram em ND (os períodos das
        pontas ficam parciais e os de fora do ciclo são descartados).
    ciclo : int, opcional
        Duração do ciclo, em dias, contando o dia da semeadura.
    ano : int, opcional
        Ano do primeiro período quando o df tem 'dia' e 'mes' sem 'ano'.
        Padrão 2023.
    a, b : float, opcional
        Coeficientes de Angström-Prescott para n/N. Padrão a = 0,29 cos(lat)
        e b = 0,52.
    T_limiar_CR, CR_frio, CR_quente : float, opcional
        Parâmetros da correção de respiração (ver `CR`). Padrão 20,0, 0,6 e 0,5.

    Retorna
    -------
    resultado : pandas.DataFrame
        Cópia do df (só os períodos do ciclo), acrescido de: 'ND', 'Qo',
        'nN', 'cTn', 'cTc', 'PPBc', 'PPBn', 'PPBp' (kg MS/ha dia), 'CR',
        'PPf_periodo' e 'PPf_acum' (kg/ha, na umidade U). A última linha
        de 'PPf_acum' é a PPf. `resultado.attrs` guarda 'PPf' e 'CIAF'.

    Exemplos
    --------
    1º decêndio de janeiro em Santa Helena-PR, milho (C4):

    >>> df = pd.DataFrame({'dia': [1], 'mes': [1], 'Tmed': [26.5], 'Qg': [22.5]})
    >>> pp = produtividade_potencial(df, lat=-24.86, rota='C4', IAF=5.0, Cc=0.35, U=13)
    >>> round(float(pp['PPBp'].iloc[0]), 0), round(pp.attrs['PPf'], 0)
    (718.0, 722.0)
    """
    if intervalo not in ('d', 'dec', 'M'):
        raise ValueError("intervalo deve ser 'd', 'dec' ou 'M'")
    if (semeadura is None) != (ciclo is None):
        raise ValueError("Informe semeadura e ciclo juntos.")
    rota = _rota(rota)

    res = df.copy()
    inicio = _inicio_periodos(res, intervalo, ano)

    # ---------- Dias de cada período no campo (ND) ----------
    if semeadura is not None:
        if inicio is None:
            raise ValueError("Para usar semeadura, o df precisa de datas ('data', DatetimeIndex ou 'dia' e 'mes').")
        sem = pd.Timestamp(semeadura).normalize()
        fim_ciclo = sem + pd.Timedelta(days=int(ciclo) - 1)
        nd_cal = (res['ND'].to_numpy(dtype=float) if 'ND' in res.columns
                  else np.array([_dias_periodo(d, intervalo) for d in inicio], dtype=float))
        fim = inicio + pd.to_timedelta(nd_cal - 1, unit='D')
        ini_campo = inicio.where(inicio > sem, sem)
        fim_campo = fim.where(fim < fim_ciclo, fim_ciclo)
        nd = np.maximum(np.asarray((fim_campo - ini_campo).days) + 1, 0)
        no_ciclo = nd > 0
        if inicio[0] > sem or fim[-1] < fim_ciclo:
            raise ValueError("O df não cobre todo o ciclo (da semeadura a semeadura + ciclo - 1).")
        res = res.loc[no_ciclo].copy()
        res['ND'] = nd[no_ciclo]
        inicio = ini_campo[no_ciclo]
    elif 'ND' not in res.columns:
        if intervalo == 'd':
            res['ND'] = 1
        elif inicio is None:
            raise ValueError("Para intervalo 'dec' ou 'M', o df precisa de datas ou da coluna 'ND'.")
        else:
            res['ND'] = [_dias_periodo(d, intervalo) for d in inicio]

    # ---------- Qo: média diária nos dias do período no campo ----------
    if 'Qo' not in res.columns:
        if inicio is None:
            raise ValueError("Sem a coluna 'Qo', o df precisa de datas para calculá-la.")
        qo = []
        for d0, nd in zip(inicio, res['ND'].astype(int)):
            nda = np.array([(d0 + pd.Timedelta(days=k)).dayofyear for k in range(nd)])
            dec = declinacao_solar(nda)
            Hn = angulo_horario_nascer(lat, dec)
            qo.append(irradiancia_extraterrestre(lat, dec, Hn, fator_correcao_distancia(nda)).mean())
        res['Qo'] = qo

    # ---------- Razão de insolação ----------
    if 'nN' in res.columns:
        res['nN'] = res['nN'].clip(0, 1)
    elif {'n', 'N'} <= set(res.columns):
        res['nN'] = (res['n'] / res['N']).clip(0, 1)
    elif 'Qg' in res.columns:
        res['nN'] = razao_insolacao(res['Qg'], res['Qo'], lat, a=a, b=b)
    else:
        raise ValueError("O df precisa de 'nN', de 'n' e 'N' ou de 'Qg'.")

    # ---------- PPBp e PPf ----------
    res['cTn'] = cTn(res['Tmed'], rota)
    res['cTc'] = cTc(res['Tmed'], rota)
    res['PPBc'] = (107.2 + 8.604 * res['Qo']) * res['cTc'] * res['nN']
    res['PPBn'] = (31.7 + 5.234 * res['Qo']) * res['cTn'] * (1 - res['nN'])
    res['PPBp'] = res['PPBc'] + res['PPBn']
    res['CR'] = CR(res['Tmed'], T_limiar=T_limiar_CR, CR_frio=CR_frio, CR_quente=CR_quente)

    ciaf = CIAF(IAF)
    res['PPf_periodo'] = res['PPBp'] * res['ND'] * ciaf * res['CR'] * Cc / (1 - 0.01 * U)
    res['PPf_acum'] = res['PPf_periodo'].cumsum()

    res.attrs['PPf'] = float(res['PPf_acum'].iloc[-1]) if len(res) else 0.0
    res.attrs['CIAF'] = ciaf
    return res


def produtividade_atingivel(PPf, ETR, ETc, ky, fase=None, metodo='produtorio'):
    """
    Produtividade atingível (PA): a PPf penalizada pelo déficit hídrico,
    pela relação de Doorenbos & Kassam (1979):

    1 - PA/PPf = ky (1 - ETR/ETc)

    - 'produtorio' (Battisti et al., 2013): aplicada fase a fase, em
      sequência; a PA ao fim de uma fase é a "PPf" da seguinte.
      PA = PPf x prod[1 - ky_i (1 - ETR_i/ETc_i)]
    - 'etapa_unica': ky do ciclo e somas de ETR e ETc no ciclo.
      PA = PPf x [1 - ky_ciclo (1 - soma ETR / soma ETc)]

    Cada fator é limitado a zero (PA nunca negativa). Fases com ETc = 0
    não penalizam.

    Parâmetros
    ----------
    PPf : float
        Produtividade potencial, em kg/ha (ver `produtividade_potencial`).
    ETR, ETc : array_like
        Evapotranspiração real e da cultura (mm), por fase ou por período
        (ex.: as colunas 'ETR' e 'ETc' de `balanco_hidrico_cultura`).
    ky : float, sequência ou dict
        Fator de resposta da produtividade ao déficit hídrico. No
        produtório, um ky por fase: sequência na ordem das fases ou dict
        {fase: ky}. Na etapa única, o ky do ciclo (escalar).
    fase : array_like, opcional
        Fase de cada período (mesmo tamanho de ETR e ETc). ETR e ETc são
        somados por fase, na ordem em que as fases aparecem. Se None, cada
        elemento de ETR e ETc já é uma fase.
    metodo : str, opcional
        'produtorio' (padrão) ou 'etapa_unica'.

    Retorna
    -------
    resultado : pandas.DataFrame
        Uma linha por fase (ou uma linha 'ciclo' na etapa única):
        'ETc', 'ETR', 'ETR/ETc', 'ky', 'fator', 'PA_inicial', 'PA_final',
        'perda' (kg/ha) e 'EC' (PA_final / PPf, acumulada). A última linha
        de 'PA_final' é a PA. `resultado.attrs` guarda 'PPf', 'PA' e 'EC'.

    Exemplos
    --------
    Milho com déficit na floração:

    >>> etc, etr = [60, 110, 120, 40], [57, 77, 102, 36]
    >>> pa = produtividade_atingivel(1000, etr, etc, ky=[0.4, 1.5, 0.5, 0.2])
    >>> round(pa.attrs['EC'], 3)
    0.489
    >>> pa1 = produtividade_atingivel(1000, etr, etc, ky=1.25, metodo='etapa_unica')
    >>> round(pa1.attrs['EC'], 3)
    0.78
    """
    ETR = np.asarray(ETR, dtype=float)
    ETc = np.asarray(ETc, dtype=float)
    if ETR.shape != ETc.shape:
        raise ValueError("ETR e ETc devem ter o mesmo tamanho.")

    if metodo == 'etapa_unica':
        if np.ndim(ky) != 0:
            raise ValueError("Na etapa única, ky deve ser o ky do ciclo (escalar).")
        tab = pd.DataFrame({'ETc': [ETc.sum()], 'ETR': [ETR.sum()]}, index=pd.Index(['ciclo'], name='fase'))
        ky_fases = [float(ky)]
    elif metodo == 'produtorio':
        if fase is None:
            tab = pd.DataFrame({'ETc': ETc, 'ETR': ETR}, index=pd.RangeIndex(len(ETc), name='fase'))
        else:
            fase = np.asarray(fase)
            if fase.shape != ETR.shape:
                raise ValueError("fase deve ter o mesmo tamanho de ETR e ETc.")
            tab = (pd.DataFrame({'fase': fase, 'ETc': ETc, 'ETR': ETR})
                   .groupby('fase', sort=False)[['ETc', 'ETR']].sum())
        if isinstance(ky, dict):
            ky_fases = [float(ky[f]) for f in tab.index]
        else:
            ky_fases = list(np.atleast_1d(np.asarray(ky, dtype=float)))
            if len(ky_fases) != len(tab):
                raise ValueError(f"ky tem {len(ky_fases)} valores para {len(tab)} fases.")
    else:
        raise ValueError("metodo deve ser 'produtorio' ou 'etapa_unica'")

    with np.errstate(divide='ignore', invalid='ignore'):
        razao = np.where(tab['ETc'] > 0, tab['ETR'] / tab['ETc'], 1.0)
    tab['ETR/ETc'] = razao
    tab['ky'] = ky_fases
    tab['fator'] = np.clip(1 - tab['ky'] * (1 - tab['ETR/ETc']), 0.0, None)

    pa = PPf * np.cumprod(tab['fator'].to_numpy())
    tab['PA_inicial'] = np.concatenate([[PPf], pa[:-1]])
    tab['PA_final'] = pa
    tab['perda'] = tab['PA_inicial'] - tab['PA_final']
    tab['EC'] = pa / PPf if PPf else np.nan

    tab.attrs['PPf'] = float(PPf)
    tab.attrs['PA'] = float(pa[-1])
    tab.attrs['EC'] = float(tab['EC'].iloc[-1])
    return tab
