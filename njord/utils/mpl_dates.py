"""
Matplotlib provides sophisticated date plotting capabilities, standing on the
shoulders of python :mod:`datetime`, the add-on modules :mod:`pytz` and
:mod:`dateutil`.


.. _date-format:

Matplotlib date format
----------------------
Matplotlib represents dates using floating point numbers specifying the number
of days since 0001-01-01 UTC, plus 1.  For example, 0001-01-01, 06:00 is 1.25,
not 0.25. Values < 1, i.e. dates before 0001-01-01 UTC are not supported.

There are a number of helper functions to convert between :mod:`datetime`
objects and Matplotlib dates:

.. currentmodule:: matplotlib.dates

.. autosummary::
   :nosignatures:

   date2num
   num2date
   num2timedelta
   epoch2num
   num2epoch
   mx2num
   drange

.. note::

   Like Python's datetime, mpl uses the Gregorian calendar for all
   conversions between dates and floating point numbers. This practice
   is not universal, and calendar differences can cause confusing
   differences between what Python and mpl give as the number of days
   since 0001-01-01 and what other software and databases yield.  For
   example, the US Naval Observatory uses a calendar that switches
   from Julian to Gregorian in October, 1582.  Hence, using their
   calculator, the number of days between 0001-01-01 and 2006-04-01 is
   732403, whereas using the Gregorian calendar via the datetime
   module we find::

     In [1]: date(2006, 4, 1).toordinal() - date(1, 1, 1).toordinal()
     Out[1]: 732401

All the Matplotlib date converters, tickers and formatters are timezone aware.
If no explicit timezone is provided, the rcParam ``timezone`` is assumend.  If
you want to use a custom time zone, pass a :class:`pytz.timezone` instance
with the tz keyword argument to :func:`num2date`, :func:`.plot_date`, and any
custom date tickers or locators you create.
See `pytz <http://pythonhosted.org/pytz/>`_ for information on :mod:`pytz` and
timezone handling.

A wide range of specific and general purpose date tick locators and
formatters are provided in this module.  See
:mod:`matplotlib.ticker` for general information on tick locators
and formatters.  These are described below.


The `dateutil module <https://dateutil.readthedocs.io/en/stable/>`_ provides
additional code to handle date ticking, making it easy to place ticks
on any kinds of dates.  See examples below.

"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import six
from six.moves import zip
import re
import time
import math
import datetime
import functools

import warnings
import logging

from dateutil.rrule import (rrule, MO, TU, WE, TH, FR, SA, SU, YEARLY,
                            MONTHLY, WEEKLY, DAILY, HOURLY, MINUTELY,
                            SECONDLY)
from dateutil.relativedelta import relativedelta
import dateutil.parser
import logging
import numpy as np

from . import units, cbook

_log = logging.getLogger(__name__)

__all__ = ('date2num', 'num2date', 'num2timedelta', 'drange', 'epoch2num',
           'num2epoch', 'mx2num', 'DateFormatter',
           'IndexDateFormatter', 'AutoDateFormatter', 'DateLocator',
           'RRuleLocator', 'AutoDateLocator', 'YearLocator',
           'MonthLocator', 'WeekdayLocator',
           'DayLocator', 'HourLocator', 'MinuteLocator',
           'SecondLocator', 'MicrosecondLocator',
           'rrule', 'MO', 'TU', 'WE', 'TH', 'FR', 'SA', 'SU',
           'YEARLY', 'MONTHLY', 'WEEKLY', 'DAILY',
           'HOURLY', 'MINUTELY', 'SECONDLY', 'MICROSECONDLY', 'relativedelta',
           'seconds', 'minutes', 'hours', 'weeks')


_log = logging.getLogger(__name__)


# Make a simple UTC instance so we don't always have to import
# pytz.  From the python datetime library docs:

class _UTC(datetime.tzinfo):
    """UTC"""

    def utcoffset(self, dt):
        return datetime.timedelta(0)

    def tzname(self, dt):
        return str("UTC")

    def dst(self, dt):
        return datetime.timedelta(0)


UTC = _UTC()


"""
Time-related constants.
"""
EPOCH_OFFSET = float(datetime.datetime(1970, 1, 1).toordinal())
JULIAN_OFFSET = 1721424.5                         # Julian date at 0001-01-01
MICROSECONDLY = SECONDLY + 1
HOURS_PER_DAY = 24.
MIN_PER_HOUR = 60.
SEC_PER_MIN = 60.
MONTHS_PER_YEAR = 12.

DAYS_PER_WEEK = 7.
DAYS_PER_MONTH = 30.
DAYS_PER_YEAR = 365.0

MINUTES_PER_DAY = MIN_PER_HOUR * HOURS_PER_DAY

SEC_PER_HOUR = SEC_PER_MIN * MIN_PER_HOUR
SEC_PER_DAY = SEC_PER_HOUR * HOURS_PER_DAY
SEC_PER_WEEK = SEC_PER_DAY * DAYS_PER_WEEK

MUSECONDS_PER_DAY = 1e6 * SEC_PER_DAY

MONDAY, TUESDAY, WEDNESDAY, THURSDAY, FRIDAY, SATURDAY, SUNDAY = (
    MO, TU, WE, TH, FR, SA, SU)
WEEKDAYS = (MONDAY, TUESDAY, WEDNESDAY, THURSDAY, FRIDAY, SATURDAY, SUNDAY)


def _to_ordinalf(dt):
    """
    Convert :mod:`datetime` or :mod:`date` to the Gregorian date as UTC float
    days, preserving hours, minutes, seconds and microseconds.  Return value
    is a :func:`float`.
    """
    # Convert to UTC
    tzi = getattr(dt, 'tzinfo', None)
    if tzi is not None:
        dt = dt.astimezone(UTC)
        tzi = UTC

    base = float(dt.toordinal())

    # If it's sufficiently datetime-like, it will have a `date()` method
    cdate = getattr(dt, 'date', lambda: None)()
    if cdate is not None:
        # Get a datetime object at midnight UTC
        midnight_time = datetime.time(0, tzinfo=tzi)

        rdt = datetime.datetime.combine(cdate, midnight_time)

        # Append the seconds as a fraction of a day
        base += (dt - rdt).total_seconds() / SEC_PER_DAY

    return base


# a version of _to_ordinalf that can operate on numpy arrays
_to_ordinalf_np_vectorized = np.vectorize(_to_ordinalf)


def _dt64_to_ordinalf(d):
    """
    Convert `numpy.datetime64` or an ndarray of those types to Gregorian
    date as UTC float.  Roundoff is via float64 precision.  Practically:
    microseconds for dates between 290301 BC, 294241 AD, milliseconds for
    larger dates (see `numpy.datetime64`).  Nanoseconds aren't possible
    because we do times compared to ``0001-01-01T00:00:00`` (plus one day).
    """

    # the "extra" ensures that we at least allow the dynamic range out to
    # seconds.  That should get out to +/-2e11 years.
    extra = d - d.astype('datetime64[s]')
    extra = extra.astype('timedelta64[ns]')
    t0 = np.datetime64('0001-01-01T00:00:00').astype('datetime64[s]')
    dt = (d.astype('datetime64[s]') - t0).astype(np.float64)
    dt += extra.astype(np.float64) / 1.0e9
    dt = dt / SEC_PER_DAY + 1.0

    NaT_int = np.datetime64('NaT').astype(np.int64)
    d_int = d.astype(np.int64)
    try:
        dt[d_int == NaT_int] = np.nan
    except TypeError:
        if d_int == NaT_int:
            dt = np.nan
    return dt


def _from_ordinalf(x, tz=None):
    """
    Convert Gregorian float of the date, preserving hours, minutes,
    seconds and microseconds.  Return value is a `.datetime`.

    The input date *x* is a float in ordinal days at UTC, and the output will
    be the specified `.datetime` object corresponding to that time in
    timezone *tz*, or if *tz* is ``None``, in the timezone specified in
    :rc:`timezone`.
    """
    if tz is None:
        tz = UTC

    ix, remainder = divmod(x, 1)
    ix = int(ix)
    if ix < 1:
        raise ValueError('Cannot convert {} to a date.  This often happens if '
                         'non-datetime values are passed to an axis that '
                         'expects datetime objects.'.format(ix))
    dt = datetime.datetime.fromordinal(ix).replace(tzinfo=UTC)

    # Since the input date `x` float is unable to preserve microsecond
    # precision of time representation in non-antique years, the
    # resulting datetime is rounded to the nearest multiple of
    # `musec_prec`. A value of 20 is appropriate for current dates.
    musec_prec = 20
    remainder_musec = int(round(remainder * MUSECONDS_PER_DAY / musec_prec)
                          * musec_prec)

    # For people trying to plot with full microsecond precision, enable
    # an early-year workaround
    if x < 30 * 365:
        remainder_musec = int(round(remainder * MUSECONDS_PER_DAY))

    # add hours, minutes, seconds, microseconds
    dt += datetime.timedelta(microseconds=remainder_musec)

    return dt.astimezone(tz)


# a version of _from_ordinalf that can operate on numpy arrays
_from_ordinalf_np_vectorized = np.vectorize(_from_ordinalf)


class strpdate2num(object):
    """
    Use this class to parse date strings to matplotlib datenums when
    you know the date format string of the date you are parsing.
    """
    def __init__(self, fmt):
        """ fmt: any valid strptime format is supported """
        self.fmt = fmt

    def __call__(self, s):
        """s : string to be converted
           return value: a date2num float
        """
        return date2num(datetime.datetime(*time.strptime(s, self.fmt)[:6]))


class bytespdate2num(strpdate2num):
    """
    Use this class to parse date strings to matplotlib datenums when
    you know the date format string of the date you are parsing.  See
    :file:`examples/misc/load_converter.py`.
    """
    def __init__(self, fmt, encoding='utf-8'):
        """
        Args:
            fmt: any valid strptime format is supported
            encoding: encoding to use on byte input (default: 'utf-8')
        """
        super(bytespdate2num, self).__init__(fmt)
        self.encoding = encoding

    def __call__(self, b):
        """
        Args:
            b: byte input to be converted
        Returns:
            A date2num float
        """
        s = b.decode(self.encoding)
        return super(bytespdate2num, self).__call__(s)


# a version of dateutil.parser.parse that can operate on nump0y arrays
_dateutil_parser_parse_np_vectorized = np.vectorize(dateutil.parser.parse)


def datestr2num(d, default=None):
    """
    Convert a date string to a datenum using
    :func:`dateutil.parser.parse`.

    Parameters
    ----------
    d : string or sequence of strings
        The dates to convert.

    default : datetime instance, optional
        The default date to use when fields are missing in *d*.
    """
    if isinstance(d, six.string_types):
        dt = dateutil.parser.parse(d, default=default)
        return date2num(dt)
    else:
        if default is not None:
            d = [dateutil.parser.parse(s, default=default) for s in d]
        d = np.asarray(d)
        if not d.size:
            return d
        return date2num(_dateutil_parser_parse_np_vectorized(d))


def date2num(d):
    """
    Convert datetime objects to Matplotlib dates.

    Parameters
    ----------
    d : `datetime.datetime` or `numpy.datetime64` or sequences of these

    Returns
    -------
    float or sequence of floats
        Number of days (fraction part represents hours, minutes, seconds, ms)
        since 0001-01-01 00:00:00 UTC, plus one.

    Notes
    -----
    The addition of one here is a historical artifact. Also, note that the
    Gregorian calendar is assumed; this is not universal practice.
    For details see the module docstring.
    """

    if hasattr(d, "values"):
        # this unpacks pandas series or dataframes...
        d = d.values

    if ((isinstance(d, np.ndarray) and np.issubdtype(d.dtype, np.datetime64))
            or isinstance(d, np.datetime64)):
        return _dt64_to_ordinalf(d)
    if not cbook.iterable(d):
        return _to_ordinalf(d)
    else:
        d = np.asarray(d)
        if not d.size:
            return d
        return _to_ordinalf_np_vectorized(d)


def julian2num(j):
    """
    Convert a Julian date (or sequence) to a Matplotlib date (or sequence).

    Parameters
    ----------
    j : float or sequence of floats
        Julian date(s)

    Returns
    -------
    float or sequence of floats
        Matplotlib date(s)
    """
    if cbook.iterable(j):
        j = np.asarray(j)
    return j - JULIAN_OFFSET


def num2julian(n):
    """
    Convert a Matplotlib date (or sequence) to a Julian date (or sequence).

    Parameters
    ----------
    n : float or sequence of floats
        Matplotlib date(s)

    Returns
    -------
    float or sequence of floats
        Julian date(s)
    """
    if cbook.iterable(n):
        n = np.asarray(n)
    return n + JULIAN_OFFSET


def num2date(x, tz=None):
    """
    Convert Matplotlib dates to `~datetime.datetime` objects.

    Parameters
    ----------
    x : float or sequence of floats
        Number of days (fraction part represents hours, minutes, seconds)
        since 0001-01-01 00:00:00 UTC, plus one.
    tz : string, optional
        Timezone of *x* (defaults to rcparams ``timezone``).

    Returns
    -------
    `~datetime.datetime` or sequence of `~datetime.datetime`
        Dates are returned in timezone *tz*.

        If *x* is a sequence, a sequence of :class:`datetime` objects will
        be returned.

    Notes
    -----
    The addition of one here is a historical artifact. Also, note that the
    Gregorian calendar is assumed; this is not universal practice.
    For details, see the module docstring.
    """
    if tz is None:
        tz = UTC
    if not cbook.iterable(x):
        return _from_ordinalf(x, tz)
    else:
        x = np.asarray(x)
        if not x.size:
            return x
        return _from_ordinalf_np_vectorized(x, tz).tolist()


def _ordinalf_to_timedelta(x):
    return datetime.timedelta(days=x)


_ordinalf_to_timedelta_np_vectorized = np.vectorize(_ordinalf_to_timedelta)


def num2timedelta(x):
    """
    Convert number of days to a `~datetime.timedelta` object.

    If *x* is a sequence, a sequence of `~datetime.timedelta` objects will
    be returned.

    Parameters
    ----------
    x : float, sequence of floats
        Number of days. The fraction part represents hours, minutes, seconds.

    Returns
    -------
    `datetime.timedelta` or list[`datetime.timedelta`]

    """
    if not cbook.iterable(x):
        return _ordinalf_to_timedelta(x)
    else:
        x = np.asarray(x)
        if not x.size:
            return x
        return _ordinalf_to_timedelta_np_vectorized(x).tolist()


def drange(dstart, dend, delta):
    """
    Return a sequence of equally spaced Matplotlib dates.

    The dates start at *dstart* and reach up to, but not including *dend*.
    They are spaced by *delta*.

    Parameters
    ----------
    dstart, dend : `~datetime.datetime`
        The date limits.
    delta : `datetime.timedelta`
        Spacing of the dates.

    Returns
    -------
    drange : `numpy.array`
        A list floats representing Matplotlib dates.

    """
    f1 = date2num(dstart)
    f2 = date2num(dend)
    step = delta.total_seconds() / SEC_PER_DAY

    # calculate the difference between dend and dstart in times of delta
    num = int(np.ceil((f2 - f1) / step))

    # calculate end of the interval which will be generated
    dinterval_end = dstart + num * delta

    # ensure, that an half open interval will be generated [dstart, dend)
    if dinterval_end >= dend:
        # if the endpoint is greated than dend, just subtract one delta
        dinterval_end -= delta
        num -= 1

    f2 = date2num(dinterval_end)  # new float-endpoint
    return np.linspace(f1, f2, num + 1)


def _close_to_dt(d1, d2, epsilon=5):
    """
    Assert that datetimes *d1* and *d2* are within *epsilon* microseconds.
    """
    delta = d2 - d1
    mus = abs(delta.total_seconds() * 1e6)
    assert mus < epsilon


def _close_to_num(o1, o2, epsilon=5):
    """
    Assert that float ordinals *o1* and *o2* are within *epsilon*
    microseconds.
    """
    delta = abs((o2 - o1) * MUSECONDS_PER_DAY)
    assert delta < epsilon


def epoch2num(e):
    """
    Convert an epoch or sequence of epochs to the new date format,
    that is days since 0001.
    """
    return EPOCH_OFFSET + np.asarray(e) / SEC_PER_DAY


def num2epoch(d):
    """
    Convert days since 0001 to epoch.  *d* can be a number or sequence.
    """
    return (np.asarray(d) - EPOCH_OFFSET) * SEC_PER_DAY


def mx2num(mxdates):
    """
    Convert mx :class:`datetime` instance (or sequence of mx
    instances) to the new date format.
    """
    scalar = False
    if not cbook.iterable(mxdates):
        scalar = True
        mxdates = [mxdates]
    ret = epoch2num([m.ticks() for m in mxdates])
    if scalar:
        return ret[0]
    else:
        return ret

def seconds(s):
    """
    Return seconds as days.
    """
    return s / SEC_PER_DAY


def minutes(m):
    """
    Return minutes as days.
    """
    return m / MINUTES_PER_DAY


def hours(h):
    """
    Return hours as days.
    """
    return h / HOURS_PER_DAY


def weeks(w):
    """
    Return weeks as days.
    """
    return w * DAYS_PER_WEEK


class DateConverter(units.ConversionInterface):
    """
    Converter for datetime.date and datetime.datetime data,
    or for date/time data represented as it would be converted
    by :func:`date2num`.

    The 'unit' tag for such data is None or a tzinfo instance.
    """

    @staticmethod
    def axisinfo(unit, axis):
        """
        Return the :class:`~matplotlib.units.AxisInfo` for *unit*.

        *unit* is a tzinfo instance or None.
        The *axis* argument is required but not used.
        """
        tz = unit

        majloc = AutoDateLocator(tz=tz)
        majfmt = AutoDateFormatter(majloc, tz=tz)
        datemin = datetime.date(2000, 1, 1)
        datemax = datetime.date(2010, 1, 1)

        return units.AxisInfo(majloc=majloc, majfmt=majfmt, label='',
                              default_limits=(datemin, datemax))

    @staticmethod
    def convert(value, unit, axis):
        """
        If *value* is not already a number or sequence of numbers,
        convert it with :func:`date2num`.

        The *unit* and *axis* arguments are not used.
        """
        return date2num(value)

    @staticmethod
    def default_units(x, axis):
        """
        Return the tzinfo instance of *x* or of its first element, or None
        """
        if isinstance(x, np.ndarray):
            x = x.ravel()

        try:
            x = cbook.safe_first_element(x)
        except (TypeError, StopIteration):
            pass

        try:
            return x.tzinfo
        except AttributeError:
            pass
        return None


units.registry[np.datetime64] = DateConverter()
units.registry[datetime.date] = DateConverter()
units.registry[datetime.datetime] = DateConverter()
