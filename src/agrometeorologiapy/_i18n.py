"""Portuguese API layer over the English core (``agrometeorologiapy.en``).

Every public function is implemented once, in English. The Portuguese
modules (``radiacao``, ``umidade``...) only declare, for each function, its
Portuguese name, how its parameters (and, where needed, string values and
DataFrame columns) translate, and its Portuguese docstring. `traduzir`
builds a wrapper with a Portuguese signature that renames the arguments,
calls the English function and translates the result back.
"""

import functools
import inspect

import numpy as np
import pandas as pd

__all__ = ["tradutor", "rotulos"]


def _inverter(mapa):
    return {v: k for k, v in mapa.items()}


def _renomear_colunas(obj, mapa):
    """Renames DataFrame columns (or dict keys) with `mapa`; other types pass through."""
    if not mapa:
        return obj
    if isinstance(obj, pd.DataFrame):
        out = obj.rename(columns=mapa)
        out.attrs = {mapa.get(k, k): v for k, v in obj.attrs.items()}
        return out
    if isinstance(obj, dict):
        return {mapa.get(k, k): v for k, v in obj.items()}
    return obj


def traduzir(func, nome, parametros=None, *, modulo, doc=None, valores=None,
             colunas_entrada=None, colunas_saida=None, fixos=None, pos=None):
    """
    Builds the Portuguese wrapper of the English function `func`.

    Parameters
    ----------
    func : callable
        English (core) function.
    nome : str
        Portuguese function name.
    parametros : dict, optional
        {Portuguese parameter name: English parameter name}. Parameters
        left out keep the English name (e.g. ``lat``).
    modulo : str
        ``__name__`` of the Portuguese module that exposes the wrapper (so
        that doctest, `help()` and `inspect` see it as defined there).
    doc : str, optional
        Portuguese docstring.
    valores : dict, optional
        {English parameter name: {Portuguese value: English value}} for
        string options (e.g. ``{'interval': {'dec': 'dekad'}}``). Defaults
        in the Portuguese signature are shown with the Portuguese value.
    colunas_entrada : dict, optional
        {Portuguese column: English column}, applied to DataFrame (or dict)
        arguments before the call.
    colunas_saida : dict, optional
        {English column: Portuguese column}, applied to the returned
        DataFrame (columns and ``attrs``) or dict (keys).
    fixos : dict, optional
        {English parameter name: value} always passed to `func` and hidden
        from the Portuguese signature (e.g. ``{'verbose': False}``, when the
        Portuguese layer prints its own messages).
    pos : callable, optional
        ``pos(result, kwargs_en)`` applied to the (already translated)
        result; its return value is returned. Used for what a column map
        cannot do, such as translating text values or printing messages.
    """
    parametros = parametros or {}
    valores = valores or {}
    en_para_pt = _inverter(parametros)
    assinatura_en = inspect.signature(func)

    nomes_en = set(assinatura_en.parameters)
    desconhecidos = (set(parametros.values()) | set(fixos or {})) - nomes_en
    if desconhecidos:
        raise ValueError(f"{nome}: parâmetros inexistentes em {func.__name__}: {sorted(desconhecidos)}")

    fixos = fixos or {}
    params_pt = []
    for p in assinatura_en.parameters.values():
        if p.name in fixos:
            continue
        padrao = p.default
        if padrao is not inspect.Parameter.empty and p.name in valores:
            padrao = _inverter(valores[p.name]).get(padrao, padrao)
        params_pt.append(p.replace(name=en_para_pt.get(p.name, p.name), default=padrao))
    assinatura_pt = assinatura_en.replace(parameters=params_pt)

    @functools.wraps(func, assigned=(), updated=())
    def wrapper(*args, **kwargs):
        ligados = assinatura_pt.bind(*args, **kwargs)
        kwargs_en = {}
        for nome_pt, valor in ligados.arguments.items():
            nome_en = parametros.get(nome_pt, nome_pt)
            if nome_en in valores and isinstance(valor, str):
                valor = valores[nome_en].get(valor, valor)
            if colunas_entrada:
                valor = _renomear_colunas(valor, colunas_entrada)
            kwargs_en[nome_en] = valor
        kwargs_en.update(fixos)
        resultado = func(**kwargs_en)
        if colunas_saida:
            if isinstance(resultado, tuple):
                resultado = tuple(_renomear_colunas(r, colunas_saida) for r in resultado)
            else:
                resultado = _renomear_colunas(resultado, colunas_saida)
        if pos is not None:
            resultado = pos(resultado, kwargs_en)
        return resultado

    wrapper.__name__ = wrapper.__qualname__ = nome
    wrapper.__module__ = modulo
    wrapper.__doc__ = doc
    wrapper.__signature__ = assinatura_pt
    wrapper.__en__ = func  # English core function
    wrapper.__fixos__ = dict(fixos)  # core arguments hidden from the Portuguese signature
    return wrapper


def rotulos(mapa):
    """Returns a `pos` function that translates class labels (a str or an
    array of str) with `mapa` ({English label: Portuguese label})."""
    def pos(resultado, kwargs_en):
        if isinstance(resultado, str):
            return mapa.get(resultado, resultado)
        return np.array([mapa.get(r, r) for r in np.ravel(resultado)], dtype=object).astype(str).reshape(np.shape(resultado))
    return pos


def tradutor(modulo):
    """Returns `traduzir` bound to the Portuguese module `modulo` (``__name__``)."""
    return functools.partial(traduzir, modulo=modulo)
