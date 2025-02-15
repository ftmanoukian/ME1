from decimal import Decimal, ROUND_HALF_UP
from statistics import mean, stdev

unitexp_conversion = {'p':-12,'n':-9,'u':-6,'m':-3,'':0,'k':3,'M':6,'G':9,'T':12}

def unitexp2int(unitexp : int | str):
    """
    Generalizes the unitexp from int | str to int. (supports 'p' through 'G' engineering exponents)
    """
    if isinstance(unitexp, str):
        unitexp = unitexp_conversion[unitexp]
    if unitexp % 3 != 0:
        raise ValueError(f'Not a valid unit exponent ({unitexp})')
    return unitexp

def unitexp2str(unitexp : int):
    """
    Returns the letter associated to the unitexp provided as an int
    """
    try:
        unitexpstr = next(k for k,v in unitexp_conversion.items() if v == unitexp)
    except StopIteration:
        raise ValueError(f'The provided unitexp ({unitexp}) does not match a valid unitexp letter')
    return unitexpstr

class Multimeter():
    unit = ''

    def __init__(self, name : str, digit_count : int):
        """
        Initializes a Multimeter object.
        Parameters:
         - name: a str representing the name/model of the instrument
         - digit_count: an integer value representing the total amount of digits (including fractional ones)
        """
        self.name = name
        self.digit_count = digit_count
        self.range_list = []

    def add_range(self, range : Decimal, unitexp : int | str, err : Decimal, n : int):
        """
        Adds the specs related to a certain scale.
        Parameters:
         - range: the maximum value of the scale (Decimal)
         - unitexp : the range's unit exponent   (int | str)
         - err: the percentual part of the error (Decimal)
         - n: the amount of error digits         (int)
        """
        unitexp = unitexp2int(unitexp)
        self.range_list.append(
            {'range': range,'unitexp': unitexp, 'err': err, 'n': n}
        )
        self.range_list.sort(key = lambda x: x['range'] * Decimal(10 ** x['unitexp'])) # reorder the list in case the ranges are provided out-of-order
        return self
    
    def range_exists(self, range : Decimal, unitexp : int | str) -> bool:
        """
        Returns True if range is specified for the instrument
        """
        unitexp = unitexp2int(unitexp)
        return any([range_data['range'] == range and range_data['unitexp'] == unitexp for range_data in self.range_list])
    
    def get_range_error_components(self, range : Decimal, unitexp : int | str) -> tuple[Decimal, int]:
        """
        Returns both error components for a given range, given the range was specified.
        Return values:
         - (percentual error, digit error)
        """
        unitexp = unitexp2int(unitexp)
        if not self.range_exists(range, unitexp):
            raise ValueError(f'The specified range ({range} {unitexp2str(unitexp)}{self.unit}) is not available for this instrument ({self.name})')
        range_data = list(filter(lambda x: x['range'] == range and x['unitexp'] == unitexp, self.range_list))[0]
        return (range_data['err'],range_data['n'])
    
    def get_best_range(self, value : Decimal, value_unitexp : int | str = 0) -> list[Decimal|int]:
        """
        Returns the lowest range that can contain the value
        """
        value_unitexp = unitexp2int(value_unitexp)
        valid_ranges = list(filter(lambda x: x['range'] * Decimal(10 ** x['unitexp'])>= value * (10 ** value_unitexp), self.range_list))
        if not any(valid_ranges):
            raise ValueError(f'The specified value ({value} {unitexp2str(value_unitexp)}{self.unit}) exceeds the instrument\'s limits ({max([range_data['range'] for range_data in self.range_list])} {self.unit}) ')
        return (valid_ranges[0]['range'],valid_ranges[0]['unitexp']) #given the ranges are ordered in ascending order, return the smaller range that fits the given value
    
    def get_displayed_value(self, value : Decimal, value_unitexp : int | str = 0, range : Decimal | None = None, range_unitexp : int | str = None) -> str:
        """
        Renders the provided value to a str as shown by the instrument.
        If range is not provided, it defaults to most accurate one for the value.
        If the value is higher than the provided range, an error is raised.
        """
        value_unitexp = unitexp2int(value_unitexp)
        if range is None or range_unitexp is None:
            range, range_unitexp = self.get_best_range(value, value_unitexp)
        range_unitexp = unitexp2int(range_unitexp)
        if not self.range_exists(range, range_unitexp):
            raise ValueError(f'The specified range ({range} {unitexp2str(range_unitexp)}{self.unit}) is not available for this instrument ({self.name})')
        else:
            if value * Decimal(10 ** value_unitexp) > range * Decimal(10 ** range_unitexp):
                return 'OL'
            #if value > range:
            #    raise ValueError(f'The specified value ({value}) exceeds the maximum for the specified range ({range})')

        value *= Decimal(10 ** (value_unitexp - range_unitexp))

        roundingexp = range.adjusted() - self.digit_count + 1
        rounded_value = str(value.quantize(Decimal(10) ** roundingexp, rounding = ROUND_HALF_UP))
        rounded_value = '0' * (self.digit_count - len(rounded_value) + rounded_value.count('.')) + rounded_value
        value_unit = unitexp2str(value_unitexp)
        return (rounded_value, value_unit)
    
    def get_indicated_counts(self, value : Decimal, value_unitexp : int | str = 0, range : Decimal | None = None, range_unitexp : int | str = None) -> int:
        """
        Returns the amount of counts displayed by the instrument for the specified value and range
        """
        displayed_value = self.get_displayed_value(value, value_unitexp, range, range_unitexp)[0]
        indicated_counts = Decimal(displayed_value.replace('.',''))
        return indicated_counts
    
    def get_error(self, value : Decimal, value_unitexp : int | str = 0, range : Decimal | None = None, range_unitexp : int | str = None) -> Decimal:
        """
        Returns the absolute error when measuring the specified value in the specified range.
        If a range is not specified, the best range is automatically selected.
        """
        if range is None or range_unitexp is None:
            range, range_unitexp = self.get_best_range(value, value_unitexp)
        elif not self.range_exists(range, range_unitexp):
            raise ValueError(f'The specified range ({range} {unitexp2str(range_unitexp)}{self.unit}) is not available for this instrument ({self.name})')
        displayed_value = self.get_displayed_value(value, value_unitexp, range, range_unitexp)[0]
        displayed_value = Decimal(displayed_value)
        indicated_counts = self.get_indicated_counts(value, value_unitexp, range, range_unitexp)
        err, n = self.get_range_error_components(range, range_unitexp)
        error = displayed_value * (err / 100 + n / indicated_counts)
        return error
    
def get_list_uncertainty(instrument : Multimeter, value_list : list[Decimal], value_list_unitexp : int | str, range : Decimal | None = None, range_unitexp : int | str = None) -> tuple[Decimal,Decimal,Decimal]:
    """
    Returns the A, B and C uncertainties for the value_list sequence measured by instrument.
    """
    if len(value_list) == 0:
        raise ValueError('An empty list was provided')
    value_list_unitexp = unitexp2int(value_list_unitexp)
    range_unitexp = unitexp2int(range_unitexp)
    if any([not isinstance(element, Decimal) for element in value_list]):
        raise TypeError('All values in value_list must be of type Decimal')
    list_mean = mean(value_list)
    list_stdev = stdev(value_list)

    u_a = list_stdev / Decimal(len(value_list)).sqrt()
    u_b = instrument.get_error(value = list_mean, value_unitexp = value_list_unitexp, range = range, range_unitexp = range_unitexp) / Decimal(3).sqrt()
    u_c = (u_a ** 2 + u_b ** 2).sqrt()

    return u_a, u_b, u_c
    
def str2declist(values : str, sep : str = ',') -> list[Decimal]:
    """
    Converts a list specified as a single str to a list of Decimal values
    """
    return list(map(Decimal, values.split(sep)))

class Voltmeter(Multimeter):
    unit = 'V'

class Ammeter(Multimeter):
    unit = 'A'

def print_measurement_result(value : Decimal, uncertainty : Decimal, unitexp : str, unit : str):
    roundingexp = (uncertainty.adjusted() - 1)
    quantized_v = value.quantize(Decimal(10) ** roundingexp, rounding = ROUND_HALF_UP)
    quantized_u = uncertainty.quantize(Decimal(10) ** roundingexp, rounding = ROUND_HALF_UP)
    print(f'({quantized_v} ± {quantized_u}) {unitexp2str(unitexp)}{unit}')

vi_instr = Voltmeter('UT61E', digit_count = 5)
vi_instr.add_range(range = Decimal('220'),  unitexp = 'm', err = Decimal('0.1'), n = 5)
vi_instr.add_range(range = Decimal('2.2'),  unitexp = '',  err = Decimal('0.1'), n = 2)
vi_instr.add_range(range = Decimal('22'),   unitexp = '',  err = Decimal('0.1'), n = 2)
vi_instr.add_range(range = Decimal('220'),  unitexp = '',  err = Decimal('0.1'), n = 2)
vi_instr.add_range(range = Decimal('1000'), unitexp = '',  err = Decimal('0.1'), n = 5)

vv_instr = Voltmeter('GDM8145', digit_count = 5)
vv_instr.add_range(range = Decimal('200'),  unitexp = 'm', err = Decimal('0.03'), n = 4)
vv_instr.add_range(range = Decimal('2'),    unitexp = '',  err = Decimal('0.03'), n = 4)
vv_instr.add_range(range = Decimal('20'),   unitexp = '',  err = Decimal('0.03'), n = 4)
vv_instr.add_range(range = Decimal('200'),  unitexp = '',  err = Decimal('0.03'), n = 4)
vv_instr.add_range(range = Decimal('1000'), unitexp = '',  err = Decimal('0.03'), n = 4)

vv = str2declist('99.96, 99.9, 99.85, 99.8, 99.74, 99.7, 99.67, 99.65, 99.61, 99.57')
vi = str2declist('99.98, 99.93, 99.85, 99.8, 99.76, 99.72, 99.68, 99.64, 99.59, 99.59')

"""
vvmed = mean(vv)
vimed = mean(vi)

uavv, ubvv, ucvv = get_list_uncertainty(vv_instr, vv, 'm', Decimal(200), 'm')
uavi, ubvi, ucvi = get_list_uncertainty(vi_instr, vi, 'm', Decimal(220), 'm')

e_cal = vimed - vvmed
u_ecal = (uavi ** 2 + uavv ** 2 + ubvv ** 2).sqrt()
U_ecal = 2 * u_ecal

err_vv = vi_instr.get_error(vimed, 'm', Decimal(220), 'm')

print("VERIFICACION")
print(f'V_i = {vimed}')
print(f'e_cal = V_i - V_v = {vimed} - {vvmed} = {e_cal}')
print(f'u_a(V_i) = {uavi}')
print(f'u_a(V_v) = {uavv}')
print(f'u_b(V_v) = {ubvv}')
print(f'u_ecal = {u_ecal}')
print(f'U_ecal = {U_ecal}')
print(f'err_vv = {err_vv}')

print(f'|e_cal| + U_ecal <= err_vv')
print(f'{abs(e_cal)} + {U_ecal} <= {err_vv}')
print(f'{abs(e_cal) + U_ecal} <= {err_vv}')

print("\n\n\nCALIBRACION")
"""