from PIL import Image

#read the image
im = Image.open("query_image.jpg")

#rotate image by 90 degrees
angle = 90
out = im.rotate(angle, expand=True)
out.save('rquery_rotate.jpg')