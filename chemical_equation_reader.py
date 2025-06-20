import re

def turn_equation_into_list(equation: str) -> tuple[list[str], list[str]] | None:
    '''
    将化学方程式转成 reactants/products 的两个列表。

    Args:
        equation (str): 化学方程式，例如 '2H2 + O2 === 2H2O'

    Returns:
        tuple[list[str], list[str]]: 反应物和生成物的两个列表，元素重复次数等于其前的系数。
        如果格式错误，返回 None。
    '''
    def parse_side(side: str) -> list[str] | None:
        compounds = side.split('+')
        result = []
        for compound in compounds:
            compound = compound.strip()
            match = re.match(r'^(\d*)([A-Za-z0-9()]+)$', compound)
            if match:
                coeff_str, formula = match.groups()
                coeff = int(coeff_str) if coeff_str else 1
                result.extend([formula] * coeff)
            else:
                return None
        return result

    # 允许一个或多个等号作为分隔符
    split_match = re.split(r'={2,}', equation)
    if len(split_match) != 2:
        return None  # 没有合法的等号分隔，或等号太多段

    left, right = split_match
    reactants = parse_side(left.strip())
    products = parse_side(right.strip())

    if reactants is None or products is None:
        return None

    return reactants, products
