{% extends 'markdown.tpl' %}

{% block input %}
<pre><code class="language-python">{{ cell.source }}</code>
</pre>
{% endblock input %}

{% block stream %}
{% if output.name == "stdout" %}
<pre>
{{ output.text.strip() | escape }}
</pre>
{% endif %}
{% endblock stream %}

{% block data_text scoped %}
<pre>
{{ output.data['text/plain'].strip() | escape }}
</pre>
{% endblock data_text %}

