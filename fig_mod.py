from PIL import Image

img = Image.open("specific.png")

print(img.info.get("dpi"))

scale = 1

new_size = (int(img.width * scale), int(img.height * scale))

img_hd = img.resize(new_size, Image.Resampling.LANCZOS)

img_hd = img_hd.convert("RGB")

img_hd.save("specific.jpeg", dpi=(900, 900))

img = Image.open("specific.jpeg")
print(img.info.get("dpi"))