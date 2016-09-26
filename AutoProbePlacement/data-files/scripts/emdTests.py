from PIL import Image;
import emd;
im1 = Image.open('/Users/joelpp/Documents/Maitrise/code/experiments/24/probes/Probe_89.png');
im2 = Image.open('/Users/joelpp/Documents/Maitrise/code/experiments/24/probes/Probe_90.png');

# print(im1.histogram()[0:255]);

print(emd.emd(range(256), range(256), im1.histogram()[0:256], im2.histogram()[0:256]));