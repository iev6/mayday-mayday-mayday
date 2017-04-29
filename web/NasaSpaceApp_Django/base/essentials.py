from django.shortcuts import render
from django.http import JsonResponse


def successResp(inputDict={}):
    success_dict = inputDict
    success_dict['status'] = True
    return JsonResponse(success_dict,status = 200)

def errorResp(status_code=400,readable_error="",extra_data={}):
    err_dict = {'status':'false','error':{'text':'','readable':readable_error}}
    msg = ''
    c = status_code
    if c == 401:
        msg = 'Unauthorized'
    elif c == 406 :
        msg = 'Data not acceptible'
    elif c == 409:
        msg = 'Conflicting data'
    elif c== 404:
        msg = 'Object not found'
    elif c == 403:
        msg = 'Forbidden'
    else :
        c = status_code
    err_dict['error']['text'] = msg
    for key in extra_data:
        err_dict[key] = extra_data[key]
    return JsonResponse(err_dict, status = c)