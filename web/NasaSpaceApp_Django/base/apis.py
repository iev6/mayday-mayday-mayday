from essentials import successResp, errorResp
from rest_framework.parsers import JSONParser
from django.shortcuts import render
import json
from dateutil import parser
from  generatecords import generate
import os

def submit(request):
	print request.META
	inputDict = request.data
	generate(inputDict)
	return successResp()

def mainCord(request):
	base_url = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
	os.system('./'+base_url+'/CARI_7_DVD/CARI7')
	os.system('cp '+base_url+'/CARI_7_DVD/PLACES.ANS '+base_url+'/web/NasaSpaceApp_nodeJS/public/json')
	os.system('mv PLACE.ANS output.csv')
	os.system('rm PLACE.ANS')
	return successResp({'data':'data'})