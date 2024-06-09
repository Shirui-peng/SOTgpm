import cdsapi

c = cdsapi.Client()

c.retrieve(
    'satellite-sea-level-global',
    {
        'variable': 'daily',
        'year': '2019',
        'month': '11',
        'day': '23',
        'version': 'vDT2021',
        'format': 'zip',
    },
    'download.zip')
