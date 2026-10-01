import struct, sys
from PIL import Image, ImageDraw, ImageFont

def render(path, out, scale=3):
    data = open(path,'rb').read()
    # header: bounds (4 ints), frame (4 ints)
    t,sz = struct.unpack_from('<II',data,0)
    bl,bt,br,bb = struct.unpack_from('<iiii',data,8)
    W, H = br-bl+1, bb-bt+1
    img = Image.new('RGB', (W*scale, H*scale), 'white')
    d = ImageDraw.Draw(img)
    objs = {}
    cur_pen = (0,0,0); cur_pen_w = 1; cur_pen_style=0
    cur_brush = None
    text_color = (0,0,0); font_h = 12; font_name='Arial'
    pos = (0,0)
    win_org=(0,0); win_ext=(W,H); vp_org=(0,0); vp_ext=(W,H)
    def tx(x,y):
        sx = vp_ext[0]/win_ext[0] if win_ext[0] else 1
        sy = vp_ext[1]/win_ext[1] if win_ext[1] else 1
        X = (x-win_org[0])*sx + vp_org[0] - bl
        Y = (y-win_org[1])*sy + vp_org[1] - bt
        return (X*scale, Y*scale)
    def getfont(h):
        try:
            return ImageFont.truetype('/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf', max(8, int(abs(h)*scale*0.95)))
        except Exception:
            return ImageFont.load_default()
    off=0
    while off < len(data):
        t,sz = struct.unpack_from('<II',data,off)
        body = data[off+8:off+sz]
        if t==9:   win_ext = struct.unpack_from('<ii',body)
        elif t==10: win_org = struct.unpack_from('<ii',body)
        elif t==11: vp_ext = struct.unpack_from('<ii',body)
        elif t==12: vp_org = struct.unpack_from('<ii',body)
        elif t==24: text_color = struct.unpack_from('<BBB',body)[:3]
        elif t==95: # EXTCREATEPEN
            ih = struct.unpack_from('<I',body)[0]
            # offBmi, cbBmi, offBits, cbBits, then LOGPENEX: penStyle, width, brushStyle, color(BBBx), hatch, numEntries
            style,width,bstyle = struct.unpack_from('<III',body,20)
            r,g,b = struct.unpack_from('<BBB',body,32)
            objs[ih] = ('pen',(r,g,b),width,style)
        elif t==38: # CREATEPEN
            ih,style,wx,wy = struct.unpack_from('<IIii',body)
            r,g,b = struct.unpack_from('<BBB',body,16)
            objs[ih]=('pen',(r,g,b),wx,style)
        elif t==39: # CREATEBRUSHINDIRECT
            ih,style = struct.unpack_from('<II',body)
            r,g,b = struct.unpack_from('<BBB',body,8)
            objs[ih]=('brush',(r,g,b),style)
        elif t==82: # EXTCREATEFONTINDIRECTW
            ih = struct.unpack_from('<I',body)[0]
            h = struct.unpack_from('<i',body,4)[0]
            name = body[4+28:4+28+64].decode('utf-16le',errors='ignore').split('\x00')[0]
            objs[ih]=('font',h,name)
        elif t==37: # SELECTOBJECT
            ih = struct.unpack_from('<I',body)[0]
            o = objs.get(ih)
            if o:
                if o[0]=='pen': cur_pen, cur_pen_w, cur_pen_style = o[1], o[2], o[3]
                elif o[0]=='brush': cur_brush = o
                elif o[0]=='font': font_h, font_name = o[1], o[2]
            else:
                # stock objects
                if ih & 0x80000000:
                    s = ih & 0xff
                    if s in (0,1,2,3,4): cur_brush=('brush',(255,255,255),0) if s==0 else ('brush',(0,0,0),0)
                    if s==5: cur_brush=None
                    if s==6: cur_pen=(255,255,255); cur_pen_w=1
                    if s==7: cur_pen=(0,0,0); cur_pen_w=1
                    if s==8: cur_pen=None
        elif t==27: pos = struct.unpack_from('<ii',body)
        elif t==54:
            p2 = struct.unpack_from('<ii',body)
            if cur_pen and cur_pen_style!=5:
                d.line([tx(*pos),tx(*p2)], fill=cur_pen, width=max(1,int(cur_pen_w*scale)))
            pos = p2
        elif t==87: # POLYLINE16
            n = struct.unpack_from('<I',body,16)[0]
            pts = struct.unpack_from('<%dh'%(2*n),body,20)
            P = [tx(pts[i],pts[i+1]) for i in range(0,2*n,2)]
            if cur_pen and len(P)>1:
                d.line(P, fill=cur_pen, width=max(1,int(cur_pen_w*scale)), joint='curve')
        elif t==43: # RECTANGLE
            l,tp,r,b2 = struct.unpack_from('<iiii',body)
            if cur_brush: d.rectangle([tx(l,tp),tx(r,b2)], fill=cur_brush[1], outline=cur_pen)
        elif t==84: # EXTTEXTOUTW
            x,y = struct.unpack_from('<ii',body,28)
            n,offs = struct.unpack_from('<II',body,36)
            s = data[off+offs:off+offs+2*n].decode('utf-16le',errors='replace')
            f = getfont(font_h)
            X,Y = tx(x,y)
            # Malvern charts: y-axis label is rotated (escapement) - check font name/escapement
            # escapement stored in LOGFONT at body offset 4+8
            esc = struct.unpack_from('<i',body,0)[0]  # not used
            d.text((X,Y), s, fill=text_color, font=f)
        if sz==0: break
        off+=sz
    img.save(out)
    return out

if __name__=='__main__':
    for i in range(1,7):
        render(f'emf/image{i}.emf', f'emf_png/dls_{i}.png')
        print('ok', i)
